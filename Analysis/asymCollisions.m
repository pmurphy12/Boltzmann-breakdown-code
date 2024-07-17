function [asymmetry] = asymCollisions(identificationString)
% Calculated number of collisions per time step for population and for
% individual cells. Additionally calculates mean free path over set
% interval and its standard deviation.

% Takes in data from 1D projection of agent-based simulation data.

    Title = 'full sim';
    SAVEFIG = true;
    SAVEDATA = true;
    
    dx = 10;

    Files = dir(fullfile('./',['*' identificationString '*.mat']));
    
    % extract information for first file to get params and length of data
    listOfVariables = who('-file',Files(1).name);
    if length(listOfVariables) ~= 3
        error(['incorrect number of variables in ' Files(1).name])
    end
    S = load(Files(1).name,listOfVariables{:});
    for ii = 1:length(listOfVariables)
        shortList1{ii} = listOfVariables{ii}(1:4);
        shortList2{ii} = listOfVariables{ii}(1:min(6,length(listOfVariables{ii})));
        shortList3{ii} = listOfVariables{ii}(1:min(10,length(listOfVariables{ii})));
    end

    testString1 = strcmp('data',shortList1);
    testString2 = strcmp('params',shortList2);
    testString3 = strcmp('collisions',shortList3);

    params = S.(listOfVariables{testString2});
    c = S.(listOfVariables{testString3});
    d = S.(listOfVariables{testString1});
    
    % Convert time to sampling indicies
    Sample = min(100,params.Nt); %How many samples taken over the simulation time.
    SampleRate = params.tend/Sample;
    t = (1:params.Nt)*params.dt;
    SampleTimes = (1:Sample)*SampleRate;
    I = zeros(size(SampleTimes));
    
    for k = 1:length(SampleTimes)
        I(k) = find(SampleTimes(k) <= t,1,'first');
    end
    
    
    collisions = cell(length(Files),1);
    collisions{1} = c(:,I)';
    
    pos = cell(length(Files),Sample);
    ori = cell(length(Files),Sample);
    temp = [d{2:end}]; % don't include initial positions
    pos(1,:) = {temp.pos};
    ori(1,:) = {temp.o};
    
    % Define some intervals based on loaded params structure.
    x = 0:dx:params.lengthGx;
    angBins = [0, pi/2, pi];
    
    
    % Proceed with the rest of the files now that params structure has been
    % loaded
    for ii = 2:length(Files)
        listOfVariables = who('-file',Files(ii).name);
        if length(listOfVariables) ~= 3
            error(['incorrect number of variables in ' Files(ii).name])
        end
        S = load(Files(ii).name,listOfVariables{:});
        for iii = 1:length(listOfVariables)
            shortList1{iii} = listOfVariables{iii}(1:4);
            shortList2{iii} = listOfVariables{iii}(1:min(6,length(listOfVariables{iii})));
            shortList3{iii} = listOfVariables{iii}(1:min(10,length(listOfVariables{iii})));
        end
        
        testString1 = strcmp('data',shortList1);
        testString3 = strcmp('collisions',shortList3);

        c = S.(listOfVariables{testString3});
        d = S.(listOfVariables{testString1});

        collisions{ii} = c(:,I)';
        temp = [d{2:end}]; % don't include initial positions
        pos(ii,:) = {temp.pos};
        ori(ii,:) = {temp.o};
    
    end
    
    collisionsM = [collisions{:}];
    posM = cell(Sample,1);
    oriM = cell(Sample,1);
    for i = 1:Sample
       posM{i} = cat(1,pos{:,i});
       oriM{i} = cat(1,ori{:,i});
    end
    
    asymmetry = zeros(size(collisionsM,1),length(x)-1);
    
    for i = 1:size(collisionsM,1)
        p = posM{i};
        Ix = discretize(p(:,1),x);
        IcL = oriM{i} >= angBins(2) & oriM{i} < angBins(3);
        IcR = oriM{i} >= angBins(1) & oriM{i} < angBins(2);
        for ii = 1:length(x)-1
            Ixi = Ix == ii;
            asymmetry(i,ii) = (sum(collisionsM(i,IcL & Ixi))-sum(collisionsM(i,IcR & Ixi)))...
                                /sum(collisionsM(i,Ixi));
        end
        
    end
    
    
    % Stats calcs
    Stats.mean_in_space = mean(asymmetry,1);
    Stats.mean_in_time = mean(asymmetry,2);
    
    %plotting
    clear h
    h = figure;
    tl = tiledlayout(1,1);
    ax1 = axes(tl);
    plot(ax1,t(I),Stats.mean_in_time)
    ax1.YLim = [-0.2 0.2];
    xlabel('time')
    ylabel('mean asymmetry (time)')
    
    ax2 = axes(tl);
    plot(ax2,x(1:end-1),Stats.mean_in_space,'Color',[0.8 0.2 0.2])
    ax2.XAxisLocation = 'top';
    ax2.YAxisLocation = 'right';
%     ax2.YLim = [-0.2 0.2];
    xlabel('x')
    ylabel('mean asymmetry (space)')
    
    ax2.Color = 'none';
    ax2.XColor = [0.8 0.2 0.2];
    ax2.YColor = [0.8 0.2 0.2];
    ax1.Box = 'off';
    ax2.Box = 'off';
    
    title(['mean asymmetry in collisions for kappa=' num2str(params.K)])
%     legend('mean asym in space','mean asym in time','Location','northwest')
    if SAVEFIG
        savefig(['collision asym means'])
        saveas(h,['collision asym means' '.png'],'png')
    end
    
    for i = [1 floor(Sample/2) Sample]
        clear h
        h = figure;
        plot(x(1:end-1),asymmetry(i,:))
        xlabel('x')
        ylabel('asymmetry in collisions')
        title(['asymmetry in collisions time-step ' num2str(i)])
    %     legend('Raw Collisions','Smoothed Collisions')
        if SAVEFIG
            savefig(['collision asym time-step ' num2str(i)])
            saveas(h,['collision asym time-step ' num2str(i) '.png'],'png')
        end
    end
    
    if SAVEDATA
        save(['Statistics for asymmetry in collisions in ' Title],'Stats','asymmetry')
    end


end
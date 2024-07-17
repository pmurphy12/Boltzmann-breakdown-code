function [Stats,C,Cjoint] = correlationAnalysisBins(tStep,identificationString,varargin)

% Takes in data from 1D projection of agent-based simulation data.

%     p = inputParser;
%     addOptional(p,'LengthX', 1, @isnumeric);
%     addOptional(p,'LengthY',1,@isnumeric);
%     addOptional(p,'save_name','',@ischar);
%     addOptional(p,'simType',' Trianglewave Test',@ischar);
%     addRequired(p,'MCdata');
%     addRequired(p,'MCparams');
    
    PLOT = true;
    SAVE = false;
    N = 2^5; % divide space into N^2 square regions
    Fact = 4;
    Edges = linspace(-pi,pi,8+1); % bin edges for orientations. correlations are calculated between these bins.
    nBins = length(Edges)-1;
    
%     TITLE = ['Local Correlations N = ' num2str(N)];
    TITLE = ['Local Independence Ratio N = ' num2str(N)];
    
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
    end

    testString1 = strcmp('data',shortList1);
    testString2 = strcmp('params',shortList2);

    params = S.(listOfVariables{testString2});
    
    data = S.(listOfVariables{testString1});
    
    % create storage for number of cells with given orientations (N1,N2)

    x = linspace(0,params.lengthGx,N/Fact+1);
    y = linspace(0,params.lengthGy,N+1);
    C = zeros(length(Files),N^2/Fact,nBins);
    Cjoint = zeros(length(Files),N^2/Fact,(nBins)^2);
    
    Ang = data{tStep}.o;
    [~,~,angI] = histcounts(Ang,Edges);
    Pos = data{tStep}.pos;
    
    %Calculate sample probabilities
    for i = 1:N/Fact
        for j = 1:N
            I = (Pos(:,1) >= x(i)) & (Pos(:,1) < x(i+1)) ...
                & (Pos(:,2) >= y(j)) & (Pos(:,2) < y(j+1));
            n = sum(I);
            
            for io = 1:nBins
                Ii = angI == io;
                angi = sum(I & Ii);
                for jo = 1:nBins
                    if io == jo         
                        if  ~isempty(angi) %Added in normalization based on conditional expectation!!!!!
                            Cjoint(1,(j-1)*N/Fact+i,(jo-1)*(nBins)+io) = angi*(angi-1)/n/(n-1);
                        else
                            Cjoint(1,(j-1)*N/Fact+i,(jo-1)*(nBins)+io) = NaN;
                        end
                    else
                        Ij = angI == jo;
                        angj = sum(I & Ij);
                        if  ~isempty(angi) && ~isempty(angj)
                            Cjoint(1,(j-1)*N/Fact+i,(jo-1)*(nBins)+io) = angi*angj/n/(n-1)/2;
                        else
                            Cjoint(1,(j-1)*N/Fact+i,(jo-1)*(nBins)+io) = NaN;
                        end                      
                    end
                end
                
                % calculate single orientation probabilities
                if n ~= 0
                    C(1,(j-1)*N/Fact+i,io) = angi/n;
                else
                    C(1,(j-1)*N/Fact+i,io) = NaN;
                end
            end
        end
    end

    % Proceed with the rest of the files now that params structure has been
    % loaded
    for ii = 2:length(Files)
        ii
        listOfVariables = who('-file',Files(ii).name);
        if length(listOfVariables) ~= 3
            error(['incorrect number of variables in ' Files(ii).name])
        end
        S = load(Files(ii).name,listOfVariables{:});
        for iii = 1:length(listOfVariables)
            shortList1{iii} = listOfVariables{iii}(1:4);
%             shortList2{iii} = listOfVariables{iii}(1:min(6,length(listOfVariables{iii})));
        end
        
        testString1 = strcmp('data',shortList1);
%         testString2 = strcmp('params',shortList2);
        
        data = S.(listOfVariables{testString1});

        Ang = data{tStep}.o;
        [~,~,angI] = histcounts(Ang,Edges);
        Pos = data{tStep}.pos;
        
        for i = 1:N/Fact
            for j = 1:N
                I = (Pos(:,1) >= x(i)) & (Pos(:,1) < x(i+1)) ...
                    & (Pos(:,2) >= y(j)) & (Pos(:,2) < y(j+1));
                n = sum(I);

                for io = 1:nBins
                    Ii = angI == io;
                    angi = sum(I & Ii);
                    for jo = 1:nBins
                        if io == jo         
                            if  ~isempty(angi) %Added in normalization based on conditional expectation!!!!!
                                Cjoint(ii,(j-1)*N/Fact+i,(jo-1)*(nBins)+io) = angi*(angi-1)/n/(n-1);
                            else
                                Cjoint(ii,(j-1)*N/Fact+i,(jo-1)*(nBins)+io) = NaN;
                            end
                        else
                            Ij = angI == jo;
                            angj = sum(I & Ij);
                            if  ~isempty(angi) && ~isempty(angj)
                                Cjoint(ii,(j-1)*N/Fact+i,(jo-1)*(nBins)+io) = angi*angj/n/(n-1)/2;
                            else
                                Cjoint(ii,(j-1)*N/Fact+i,(jo-1)*(nBins)+io) = NaN;
                            end                      
                        end
                    end

                    % calculate single orientation probabilities
                    if n ~= 0
                        C(ii,(j-1)*N/Fact+i,io) = angi/n;
                    else
                        C(ii,(j-1)*N/Fact+i,io) = NaN;
                    end
                end
            end
        end
    
    end
    
    
    Stats.N1Marginal = squeeze(mean(C,1,'omitnan'));
    Stats.N1N2joint = squeeze(mean(Cjoint,1,'omitnan'));
%     Stats.m = angI(1)*Stats.N1Marginal'+angI(2)*Stats.N2Marginal';
% 
%     Stats.correlation = (angI(1)^2*Stats.N1N2joint(:,1)+angI(1)*angI(2)*Stats.N1N2joint(:,2)...
%                         +angI(1)*angI(2)*Stats.N1N2joint(:,3)+angI(2)^2*Stats.N1N2joint(:,4)...
%                         -Stats.m.^2)...
%                         ./((angI(1)-Stats.m).^2.*Stats.N1Marginal'+(angI(2)-Stats.m).^2.*Stats.N2Marginal');
%     Stats.fixedSpaceCorr = reshape(Stats.correlation,[N/Fact,N]);
    for i = 1:nBins
        for j = 1:nBins
            Stats.independence(:,:,(j-1)*nBins+i) = reshape(Stats.N1N2joint(:,(j-1)*nBins+i)-Stats.N1Marginal(:,i).*Stats.N1Marginal(:,j),[N/Fact,N]);
            Stats.independenceRatio(:,:,(j-1)*nBins+i) = reshape((Stats.N1N2joint(:,(j-1)*nBins+i)-Stats.N1Marginal(:,i).*Stats.N1Marginal(:,j))./(Stats.N1Marginal(:,i).*Stats.N1Marginal(:,j)),[N/Fact,N]);
        end
    end

%     [Stats.correlation,Stats.pvals] = corr(numC1,numC2);
%     Stats.fixedSpaceCorr = reshape(diag(Stats.correlation),[N,N]);
%     Stats.fixedSpacePvals = reshape(diag(Stats.pvals),[N,N]);

    save(['Statistical correlations and marginals for time step ' num2str(tStep) ' and spatial discretization ' num2str(N)],'Stats','C','Cjoint')
    
    if PLOT
        [X,Y] = meshgrid(x(1:end-1),y(1:end-1));
        h = figure('Position',[300 100 300+2*30 300+2*300]); 
%         surf(X,Y,Stats.fixedSpaceCorr');
        surf(X,Y,Stats.independence(:,:,3)');
        view([0 90])
        colorbar
        caxis([-1/nBins/(nBins-1) 1/nBins/(nBins-1)])
        xlabel('x')
        ylabel('y')
        title([TITLE ' Timestep ' num2str(tStep)])
        
        if SAVE
            savefig([TITLE ' 3 Timestep ' num2str(tStep)])
            saveas(h,[TITLE ' 3 Timestep ' num2str(tStep) '.png'],'png')
        end
    end

end

% TITLE = 'Local Independence 11'
% if PLOT
% [X,Y] = meshgrid(x(1:end-1),y(1:end-1));
% h = figure('Position',[300 100 300+2*30 300+2*300]);
% surf(X,Y,Stats.independence11');
% view([0 90])
% colorbar
% caxis([-1 1])
% xlabel('x')
% ylabel('y')
% title([TITLE ' Timestep ' num2str(tStep)])
% if SAVE
% savefig([TITLE ' Timestep ' num2str(tStep)])
% saveas(h,[TITLE ' Timestep ' num2str(tStep) '.png'],'png')
% end
% end
% TITLE = 'Local Independence 12'
% if PLOT
% [X,Y] = meshgrid(x(1:end-1),y(1:end-1));
% h = figure('Position',[300 100 300+2*30 300+2*300]);
% surf(X,Y,Stats.independence12');
% view([0 90])
% colorbar
% caxis([-1 1])
% xlabel('x')
% ylabel('y')
% title([TITLE ' Timestep ' num2str(tStep)])
% if SAVE
% savefig([TITLE ' Timestep ' num2str(tStep)])
% saveas(h,[TITLE ' Timestep ' num2str(tStep) '.png'],'png')
% end
% end
% TITLE = 'Local Independence 21'
% if PLOT
% [X,Y] = meshgrid(x(1:end-1),y(1:end-1));
% h = figure('Position',[300 100 300+2*30 300+2*300]);
% surf(X,Y,Stats.independence21');
% view([0 90])
% colorbar
% caxis([-1 1])
% xlabel('x')
% ylabel('y')
% title([TITLE ' Timestep ' num2str(tStep)])
% if SAVE
% savefig([TITLE ' Timestep ' num2str(tStep)])
% saveas(h,[TITLE ' Timestep ' num2str(tStep) '.png'],'png')
% end
% end
% TITLE = 'Local Independence 22'
% if PLOT
% [X,Y] = meshgrid(x(1:end-1),y(1:end-1));
% h = figure('Position',[300 100 300+2*30 300+2*300]);
% surf(X,Y,Stats.independence22');
% view([0 90])
% colorbar
% caxis([-1 1])
% xlabel('x')
% ylabel('y')
% title([TITLE ' Timestep ' num2str(tStep)])
% if SAVE
% savefig([TITLE ' Timestep ' num2str(tStep)])
% saveas(h,[TITLE ' Timestep ' num2str(tStep) '.png'],'png')
% end
% end
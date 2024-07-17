function [Stats,C1_num,C2_num,Cjoint_num,params] = correlationAnalysis_corrected_nonconditional2(tStep,identificationString,varargin)

% Takes in data from 1D projection of agent-based simulation data.

    p = inputParser;
    addOptional(p,'N', 5, @isnumeric);
    addOptional(p,'Fact', 1, @isnumeric);
    addOptional(p,'SAVEFIG', false, @islogical);
%     addOptional(p,'LengthX', 1, @isnumeric);
%     addOptional(p,'LengthY',1,@isnumeric);
%     addOptional(p,'save_name','',@ischar);
%     addOptional(p,'simType',' Trianglewave Test',@ischar);
%     addRequired(p,'MCdata');
%     addRequired(p,'MCparams');
    parse(p,varargin{:}); 
    
    PLOT = false;
    SAVE = false;
    SAVEDATA = false;
    N = p.Results.N; % divide space into N^2 square regions
    Fact = p.Results.Fact; % divide region into rectangles with width-to-heigth ration Fact
%     TITLE = ['Local Correlations N = ' num2str(N)];
    TITLE = ['Local Independence Ratio N = ' num2str(N)];
    
    Files = dir(fullfile('./',['*' identificationString '*.mat']));
    Files = Files(1:100); % comment out to do full analysis
    
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
    Ang = data{tStep}.o;
    ang = unique(Ang);
    
    % create storage for number of cells with given orientations (N1,N2)

    x = linspace(0,params.lengthGx,N/Fact+1);
    y = linspace(0,params.lengthGy,N+1);
    C1_num = zeros(length(Files),N^2/Fact);
    C2_num = zeros(length(Files),N^2/Fact);
    Cjoint_num = zeros(length(Files),N^2/Fact,4);
    C_denom = zeros(length(Files),N^2/Fact);
    C_denom_pairs = zeros(length(Files),N^2/Fact);
    
        
    %Calculate random variable pair (N1,N2)
    for ii = 1:length(Files)
%         ii
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
        Pos = data{tStep}.pos;
        
    for i = 1:N/Fact
        for j = 1:N
            I = (Pos(:,1) > x(i)) & (Pos(:,1) < x(i+1)) & (Pos(:,2) > y(j)) & (Pos(:,2) < y(j+1));
            AngTemp = Ang(I);
            ang1 = AngTemp(AngTemp == ang(1));
            ang2 = AngTemp(AngTemp == ang(2));
            n = length(AngTemp);
            A = AngTemp*AngTemp';
            A(1:n+1:n^2) = NaN;
            if n ~= 0
                C_denom(ii,(j-1)*N/Fact+i) = n;
                C_denom_pairs(ii,(j-1)*N/Fact+i) = n*(n-1);
                C1_num(ii,(j-1)*N/Fact+i) = length(ang1);%/params.N;
                C2_num(ii,(j-1)*N/Fact+i) = length(ang2);%/params.N;
                testCases = [ang(1)^2, ang(1)*ang(2), ang(2)^2];
                Cjoint_num(ii,(j-1)*N/Fact+i,1) = sum(A == testCases(1),'all');%/params.N/(params.N-1);
                Cjoint_num(ii,(j-1)*N/Fact+i,2) = sum(A == testCases(2),'all')/2;%/params.N/(params.N-1);
                Cjoint_num(ii,(j-1)*N/Fact+i,3) = sum(A == testCases(2),'all')/2;%/params.N/(params.N-1);
                Cjoint_num(ii,(j-1)*N/Fact+i,4) = sum(A == testCases(3),'all');%/params.N/(params.N-1);
                if length(AngTemp)>5
                    5+1;
                end
            else
                C_denom(ii,(j-1)*N/Fact+i) = NaN;
                C_denom_pairs(1,(j-1)*N/Fact+i) = NaN;
                C1_num(ii,(j-1)*N/Fact+i) = NaN;
                C2_num(ii,(j-1)*N/Fact+i) = NaN;
                Cjoint_num(ii,(j-1)*N/Fact+i,1) = NaN;
                Cjoint_num(ii,(j-1)*N/Fact+i,2) = NaN;
                Cjoint_num(ii,(j-1)*N/Fact+i,3) = NaN;
                Cjoint_num(ii,(j-1)*N/Fact+i,4) = NaN;
            end
            
        end
    end
    
    end
    
    %Calculate conditional probabilities
%     Cjoint_num(:,:,1) = Cjoint_num(:,:,1)./C_denom_pairs;
%     Cjoint_num(:,:,2) = Cjoint_num(:,:,2)./C_denom_pairs;
%     Cjoint_num(:,:,3) = Cjoint_num(:,:,3)./C_denom_pairs;
%     Cjoint_num(:,:,4) = Cjoint_num(:,:,4)./C_denom_pairs;
%     C1_num = C1_num./C_denom;
%     C2_num = C2_num./C_denom;
    % take averages over simulations
    Stats.N1Marginal = mean(C1_num,1,'omitnan');%./Stats.denom;
    Stats.N2Marginal = mean(C2_num,1,'omitnan');%./Stats.denom;
    Stats.N1N2joint = squeeze(mean(Cjoint_num,1,'omitnan'));
%     Stats.N1N2joint(:,1) = Stats.N1N2joint(:,1)./Stats.denom_pairs';
%     Stats.N1N2joint(:,2) = Stats.N1N2joint(:,2)./Stats.denom_pairs';
%     Stats.N1N2joint(:,3) = Stats.N1N2joint(:,3)./Stats.denom_pairs';
%     Stats.N1N2joint(:,4) = Stats.N1N2joint(:,4)./Stats.denom_pairs';

    % calculate probabilities of rod (or pair) being in subregion ij
    Stats.denom = mean(C_denom,1,'omitnan');%/params.N;
    Stats.denom_pairs = mean(C_denom_pairs,1,'omitnan');%/params.N/(params.N);
    
    % pooled statistics
    Stats.N1Marginal_pooled = sum(C1_num,1,'omitnan')./sum(C_denom,1,'omitnan');
    Stats.N2Marginal_pooled = sum(C2_num,1,'omitnan')./sum(C_denom,1,'omitnan');
    temp = squeeze(sum(Cjoint_num,1,'omitnan'));
    temp2 = sum(C_denom_pairs,1,'omitnan')';
    Stats.N1N2joint_pooled(:,1) = temp(:,1)./temp2;
    Stats.N1N2joint_pooled(:,2) = temp(:,2)./temp2;
    Stats.N1N2joint_pooled(:,3) = temp(:,3)./temp2;
    Stats.N1N2joint_pooled(:,4) = temp(:,4)./temp2;
    
    % Calculate correction factor for conditional probabilities
    Stats.num_num_pair_ratio = Stats.denom.^2./Stats.denom_pairs;
    
    % calculate other statistics
    Stats.m = ang(1)*Stats.N1Marginal'+ang(2)*Stats.N2Marginal';

    Stats.correlation = (ang(1)^2*Stats.N1N2joint(:,1)+ang(1)*ang(2)*Stats.N1N2joint(:,2)...
                        +ang(1)*ang(2)*Stats.N1N2joint(:,3)+ang(2)^2*Stats.N1N2joint(:,4)...
                        -Stats.m.^2)...
                        ./((ang(1)-Stats.m).^2.*Stats.N1Marginal'+(ang(2)-Stats.m).^2.*Stats.N2Marginal');
    Stats.fixedSpaceCorr = reshape(Stats.correlation,[N/Fact,N]);
    
    Stats.independence11 = reshape(Stats.N1N2joint(:,1)-(Stats.N1Marginal.*(Stats.N1Marginal))',[N/Fact,N]);
    Stats.independence12 = reshape(Stats.N1N2joint(:,2)-(Stats.N1Marginal.*Stats.N2Marginal)',[N/Fact,N]);
    Stats.independence21 = reshape(Stats.N1N2joint(:,3)-(Stats.N1Marginal.*Stats.N2Marginal)',[N/Fact,N]);
    Stats.independence22 = reshape(Stats.N1N2joint(:,4)-(Stats.N2Marginal.*(Stats.N2Marginal))',[N/Fact,N]);
    
    Stats.independenceRatio11 = reshape((Stats.N1N2joint(:,1)-(Stats.N1Marginal.^2)')./(Stats.N1Marginal.^2'),[N/Fact,N]);
    Stats.independenceRatio12 = reshape((Stats.N1N2joint(:,2)-(Stats.N1Marginal.*Stats.N2Marginal)')./(Stats.N1Marginal'.*Stats.N2Marginal'),[N/Fact,N]);
    Stats.independenceRatio21 = reshape((Stats.N1N2joint(:,3)-(Stats.N1Marginal.*Stats.N2Marginal)')./(Stats.N1Marginal'.*Stats.N2Marginal'),[N/Fact,N]);
    Stats.independenceRatio22 = reshape((Stats.N1N2joint(:,4)-(Stats.N2Marginal.^2)')./(Stats.N2Marginal.^2'),[N/Fact,N]);

%     [Stats.correlation,Stats.pvals] = corr(numC1,numC2);
%     Stats.fixedSpaceCorr = reshape(diag(Stats.correlation),[N,N]);
%     Stats.fixedSpacePvals = reshape(diag(Stats.pvals),[N,N]);

    if SAVEDATA
        save(['test Statistical conditional correlations and marginals for time step ' num2str(tStep) ' and spatial discretization ' num2str(N)],'Stats','C1_num','C2_num','Cjoint_num','C_denom')
    end
        
    if PLOT
        [X,Y] = meshgrid(x(1:end-1),y(1:end-1));
%         h = figure('Position',[300 100 300+2*30 300+2*300]); 
        h = figure('Position',[300 100 600 600]); 
%         surf(X,Y,Stats.fixedSpaceCorr');
        surf(X,Y,Stats.independenceRatio12');
        view([0 90])
        colorbar
        caxis([-1 1])
        xlabel('x')
        ylabel('y')
        title([TITLE ' Timestep ' num2str(tStep)])
        
        if SAVE
            savefig([TITLE ' Timestep ' num2str(tStep)])
            saveas(h,[TITLE ' Timestep ' num2str(tStep) '.png'],'png')
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
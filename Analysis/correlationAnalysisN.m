function [Stats,C1,C2,Cjoint,params] = correlationAnalysisN(tStep,identificationString,varargin)

% Takes in data from 1D projection of agent-based simulation data.

    p = inputParser;
    addOptional(p,'N', 2^5, @isnumeric);
    addOptional(p,'Fact',4,@isnumeric);
    addOptional(p,'PLOT',true,@islogical);
    addOptional(p,'SAVEFIG',false,@islogical);
    addOptional(p,'SAVEDATA',false,@islogical);
%     addRequired(p,'TITLE',@ischar);
    parse(p,varargin{:}); 
    
    PLOT = p.Results.PLOT;
    SAVEFIG = p.Results.SAVEFIG;
    SAVEDATA = p.Results.SAVEDATA;
    N = p.Results.N; % divide space into N^2 square regions
    Fact = p.Results.Fact;
    TITLE = ['Local independence for discretization = ' num2str(N)];
    TITLE2 = ['Local independence ratio for discretization = ' num2str(N)];
    
    Files = dir(fullfile('./',['*' identificationString '*.mat']));
    
    % extract information for first file to get params and length of data
    listOfVariables = who('-file',Files(1).name);
    if length(listOfVariables) == 1 || length(listOfVariables) > 3
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
    C1 = zeros(length(Files),N^2/Fact);
    C2 = zeros(length(Files),N^2/Fact);
    Cjoint = zeros(length(Files),N^2/Fact,4);
    
    Ang = data{tStep}.o;
    ang = unique(Ang);
    Pos = data{tStep}.pos;
    
    %Calculate random variable pair (N1,N2)
    for i = 1:N/Fact
        for j = 1:N
            I = (Pos(:,1) >= x(i)) & (Pos(:,1) < x(i+1)) & (Pos(:,2) >= y(j)) & (Pos(:,2) < y(j+1));
            AngTemp = Ang(I);
            ang1 = AngTemp(AngTemp == ang(1));
            ang2 = AngTemp(AngTemp == ang(2));
            n = length(AngTemp);
            A = AngTemp*AngTemp';
            A(1:n+1:n^2) = NaN;
            if n ~= 0 %Added in normalization based on conditional expectation!!!!!
                C1(1,(j-1)*N/Fact+i) = length(ang1);
                C2(1,(j-1)*N/Fact+i) = length(ang2);
                testCases = [ang(1)^2, ang(1)*ang(2), ang(2)^2];
                Cjoint(1,(j-1)*N/Fact+i,1) = sum(A == testCases(1),'all');
                Cjoint(1,(j-1)*N/Fact+i,2) = sum(A == testCases(2),'all')/2;
                Cjoint(1,(j-1)*N/Fact+i,3) = sum(A == testCases(2),'all')/2;
                Cjoint(1,(j-1)*N/Fact+i,4) = sum(A == testCases(3),'all');
            else
                C1(1,(j-1)*N/Fact+i) = NaN;
                C2(1,(j-1)*N/Fact+i) = NaN;
                Cjoint(1,(j-1)*N/Fact+i,1) = NaN;
                Cjoint(1,(j-1)*N/Fact+i,2) = NaN;
                Cjoint(1,(j-1)*N/Fact+i,3) = NaN;
                Cjoint(1,(j-1)*N/Fact+i,4) = NaN;
            end
            
        end
    end

    % Proceed with the rest of the files now that params structure has been
    % loaded
    for ii = 2:length(Files)
        ii
        listOfVariables = who('-file',Files(ii).name);
        if length(listOfVariables) == 1 || length(listOfVariables) > 3
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
                C1(ii,(j-1)*N/Fact+i) = length(ang1);
                C2(ii,(j-1)*N/Fact+i) = length(ang2);
                testCases = [ang(1)^2, ang(1)*ang(2), ang(2)^2];
                Cjoint(ii,(j-1)*N/Fact+i,1) = sum(A == testCases(1),'all');
                Cjoint(ii,(j-1)*N/Fact+i,2) = sum(A == testCases(2),'all')/2;
                Cjoint(ii,(j-1)*N/Fact+i,3) = sum(A == testCases(2),'all')/2;
                Cjoint(ii,(j-1)*N/Fact+i,4) = sum(A == testCases(3),'all');
                if length(AngTemp)>5
                    5+1;
                end
            else
                C1(ii,(j-1)*N/Fact+i) = NaN;
                C2(ii,(j-1)*N/Fact+i) = NaN;
                Cjoint(ii,(j-1)*N/Fact+i,1) = NaN;
                Cjoint(ii,(j-1)*N/Fact+i,2) = NaN;
                Cjoint(ii,(j-1)*N/Fact+i,3) = NaN;
                Cjoint(ii,(j-1)*N/Fact+i,4) = NaN;
            end
            
        end
    end
    
    end
    
    
    Stats.N1Marginal = mean(C1,1,'omitnan');
    Stats.N2Marginal = mean(C2,1,'omitnan');
    Stats.N1N2joint = squeeze(mean(Cjoint,1,'omitnan'));
    Stats.m = ang(1)*Stats.N1Marginal'+ang(2)*Stats.N2Marginal';

    Stats.correlation = (ang(1)^2*Stats.N1N2joint(:,1)+ang(1)*ang(2)*Stats.N1N2joint(:,2)...
                        +ang(1)*ang(2)*Stats.N1N2joint(:,3)+ang(2)^2*Stats.N1N2joint(:,4)...
                        -Stats.m.^2)...
                        ./((ang(1)-Stats.m).^2.*Stats.N1Marginal'+(ang(2)-Stats.m).^2.*Stats.N2Marginal');
    Stats.fixedSpaceCorr = reshape(Stats.correlation,[N/Fact,N]);
    Stats.independence11 = reshape(Stats.N1N2joint(:,1)-Stats.N1Marginal.^2',[N/Fact,N]);
    Stats.independence12 = reshape(Stats.N1N2joint(:,2)-Stats.N1Marginal'.*Stats.N2Marginal',[N/Fact,N]);
    Stats.independence21 = reshape(Stats.N1N2joint(:,3)-Stats.N1Marginal'.*Stats.N2Marginal',[N/Fact,N]);
    Stats.independence22 = reshape(Stats.N1N2joint(:,4)-Stats.N2Marginal.^2',[N/Fact,N]);
    
    Stats.independenceRatio11 = reshape((Stats.N1N2joint(:,1)-Stats.N1Marginal.^2')./(Stats.N1N2joint(:,1)),[N/Fact,N]);
    Stats.independenceRatio12 = reshape((Stats.N1N2joint(:,2)-Stats.N1Marginal'.*Stats.N2Marginal')./(Stats.N1N2joint(:,2)),[N/Fact,N]);
    Stats.independenceRatio21 = reshape((Stats.N1N2joint(:,3)-Stats.N1Marginal'.*Stats.N2Marginal')./(Stats.N1N2joint(:,3)),[N/Fact,N]);
    Stats.independenceRatio22 = reshape((Stats.N1N2joint(:,4)-Stats.N2Marginal.^2')./(Stats.N1N2joint(:,4)),[N/Fact,N]);

%     [Stats.correlation,Stats.pvals] = corr(numC1,numC2);
%     Stats.fixedSpaceCorr = reshape(diag(Stats.correlation),[N,N]);
%     Stats.fixedSpacePvals = reshape(diag(Stats.pvals),[N,N]);
    
    if SAVEDATA
        save(['Statistical correlations and marginals for time step ' num2str(tStep) ' and spatial discretization ' num2str(N)],'Stats','C1','C2','Cjoint')
    end

    if PLOT
        [X,Y] = meshgrid(x(1:end-1),y(1:end-1));
%         h = figure('Position',[300 50 300+2*30 50+2*300]); 
        h = figure('Position',[300 50 600 600]); 
%         surf(X,Y,Stats.fixedSpaceCorr');
        surf(X,Y,Stats.independence12');
        view([0 90])
        colorbar
        caxis([-1 1]*params.N/(N^2/4/4))
        xlabel('x')
        ylabel('y')
        title([TITLE ' Timestep ' num2str(tStep)])
        
        if SAVEFIG
            savefig([TITLE ' timestep ' num2str(tStep)])
            saveas(h,[TITLE ' timestep ' num2str(tStep) '.png'],'png')
        end
    end
    
    
    if PLOT
        [X,Y] = meshgrid(x(1:end-1),y(1:end-1));
%         h2 = figure('Position',[300 50 300+2*30 50+2*300]);
        h2 = figure('Position',[600 50 600 600]); 
%         surf(X,Y,Stats.fixedSpaceCorr');
        surf(X,Y,Stats.independenceRatio12');
        view([0 90])
        colorbar
        caxis([-1 1])
        xlabel('x')
        ylabel('y')
        title([TITLE2 ' timestep ' num2str(tStep)])
        
        if SAVEFIG
            savefig([TITLE2 ' timestep ' num2str(tStep)])
            saveas(h2,[TITLE2 ' timestep ' num2str(tStep) '.png'],'png')
        end
    end

end

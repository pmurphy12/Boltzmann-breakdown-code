function [clusterSizes1,clusterSizes2,Stats] = clusterAnalysisAxis(identificationString,varargin)
%CLUSTERANALYSIS Take sim data and parameters and calculates clusters for
%each time step. Additionally calculates distribution of cluster sizes.


    PLOT = true;
    SAVE = true;
    SAVEFIG = true;
    TITLE = ['Clusters in final time'];
    nBars = 25;
    
    % NOTE: change params.dt below to match min distance for different time steps
    stepScale = 1;
    
    Files = dir(fullfile('./',['*' identificationString '*.mat']));
    
    minN = 2; %minimum number of neightbors to count in a cluster
    
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
    
    
    clusterSizes1 = cell(length(Files),length(data));
    clusterSizes2 = cell(length(Files),length(data));
    clusterLength1 = cell(length(Files),length(data));
    clusterLength2 = cell(length(Files),length(data));
    EDGES = -1.5:1:(params.N/minN+0.5);
    EDGEScnumb = 0.5:1:params.N+0.5;
    
    
    parfor i = 1:length(data)
        
        ang = unique(data{i}.o);
        
        for k = ang'

            X = data{i}.pos(data{i}.o == k,:);
%             idx = dbscan(X,params.len/1.9+params.dt*params.vbar*stepScale,minN);
            idx = dbscan(X,params.len*1.05,minN);
            if k == ang(1)
                clusterSizes1{1,i} = histcounts(idx,EDGES);
            elseif k == ang(2)
                clusterSizes2{1,i} = histcounts(idx,EDGES);
            end
            
            clusterIds = setdiff(unique(idx),-1);
            clusterLengths = zeros(length(clusterIds),1);
            if ~isempty(clusterIds)
                kk = 1;
                for ii = clusterIds'
                    D = pdist2(X(idx == ii,:),X(idx == ii,:));
                    clusterLengths(kk) = max(triu(D,1),[],'all')+params.len;
                    kk = kk+1;
                end
            end
            if k == ang(1)
                clusterLength1{1,i} = clusterLengths';
            elseif k == ang(2)
                clusterLength2{1,i} = clusterLengths';
            end
            
            
        end
        
    end
    
    
    % repeate for the rest of the datasets
    for f = 2:length(Files)
        f 
        listOfVariables = who('-file',Files(f).name);
        if length(listOfVariables) ~= 3
            error(['incorrect number of variables in ' Files(f).name])
        end
        S = load(Files(f).name,listOfVariables{:});
        for ii = 1:length(listOfVariables)
            shortList1{ii} = listOfVariables{ii}(1:4);
            shortList2{ii} = listOfVariables{ii}(1:min(6,length(listOfVariables{ii})));
        end

        testString1 = strcmp('data',shortList1);
        testString2 = strcmp('params',shortList2);

        params = S.(listOfVariables{testString2});

        data = S.(listOfVariables{testString1});


        parfor i = 1:length(data)

            ang = unique(data{i}.o);

            for k = ang'

                X = data{i}.pos(data{i}.o == k,:);
                
%                idx = dbscan(X,params.len/1.9+params.dt*params.vbar*stepScale,minN);
                idx = dbscan(X,params.len*1.05,minN);
                if k == ang(1)
                    clusterSizes1{f,i} = histcounts(idx,EDGES);
                elseif k == ang(2)
                    clusterSizes2{f,i} = histcounts(idx,EDGES);
                end
                                
                clusterIds = setdiff(unique(idx),-1);
                clusterLengths = zeros(length(clusterIds),1);
                kk = 1;
                if ~isempty(clusterIds)
                    for ii = clusterIds'
                        D = pdist2(X(idx == ii,:),X(idx == ii,:));
                        clusterLengths(kk) = max(triu(D,1),[],'all')+params.len;
                        kk = kk+1;
                    end
                end
%                 clusterLengths(kk) = params.len; %add single cells to length profile
    
                if k == ang(1)
                    clusterLength1{f,i} = clusterLengths';
                elseif k == ang(2)
                    clusterLength2{f,i} = clusterLengths';
                end

            end

        end
        
        if f == length(Files) && PLOT
            ang = unique(data{end}.o);

            X = data{end}.pos(data{end}.o == ang(1),:);
%             idx = dbscan(X,params.len/1.9+params.dt*params.vbar*stepScale,minN);
            idx = dbscan(X,params.len*1.05,minN);
            h = figure; gscatter(X(:,1),X(:,2),idx)
            title([TITLE ' for angle ' num2str(1)])
            legend('off')
            if SAVEFIG
                savefig(['Scatter plot of groups in final time step for sim ang 1'])
                saveas(h,['Scatter plot of groups in final time step for last sim ang 1.png'],'png')
            end

            X = data{end}.pos(data{end}.o == ang(2),:);
            idx = dbscan(X,params.len/1.9+params.dt*params.vbar*stepScale,minN);
            h2 = figure; gscatter(X(:,1),X(:,2),idx)
            title([TITLE ' for angle ' num2str(2)])
            legend('off')
            if SAVEFIG
                savefig(['Scatter plot of groups in final time step for sim ang 2'])
                saveas(h2,['Scatter plot of groups in final time step for last sim ang 2.png'],'png')
            end
        end
        
        
    end
    
    % add data from different trials
    clusterDist1 = zeros(length(EDGEScnumb)-1,length(data));
    clusterDist2 = zeros(length(EDGEScnumb)-1,length(data));
    
    for i = 1:length(data)
        for j = 1:length(Files)
            clusterDist1(1,i) = clusterDist1(1,i) + clusterSizes1{j,i}(1);
            clusterDist2(1,i) = clusterDist2(1,i) + clusterSizes2{j,i}(1);
            clusterDist1(2:end,i) = clusterDist1(2:end,i) + histcounts(clusterSizes1{j,i}(2:end),EDGEScnumb(2:end))';
            clusterDist2(2:end,i) = clusterDist2(2:end,i) + histcounts(clusterSizes2{j,i}(2:end),EDGEScnumb(2:end))';
        end
    end
    
    Stats.meanClusterSize1 = zeros(1,length(data));
    Stats.meanClusterSize2 = zeros(1,length(data));
    Stats.stdClusterSize1 = zeros(1,length(data));
    Stats.stdClusterSize2 = zeros(1,length(data));
    % statistics for clusters not counting individual cells (or pairs technically)
    Stats.meanClusterSizeGroups1 = zeros(1,length(data));
    Stats.meanClusterSizeGroups2 = zeros(1,length(data));
    Stats.stdClusterSizeGroups1 = zeros(1,length(data));
    Stats.stdClusterSizeGroups2 = zeros(1,length(data));
    
    % statistics for cluster max diameter
    Stats.meanClusterLength1 = zeros(1,length(data));
    Stats.meanClusterLength2 = zeros(1,length(data));
    Stats.stdClusterLength1 = zeros(1,length(data));
    Stats.stdClusterLength2 = zeros(1,length(data));
    for i = 1:length(data)
        Stats.meanClusterSize1(i) = sum(clusterDist1(:,i).*(1:length(EDGEScnumb)-1)')/sum(clusterDist1(:,i));
        Stats.meanClusterSize2(i) = sum(clusterDist2(:,i).*(1:length(EDGEScnumb)-1)')/sum(clusterDist2(:,i));
        Stats.stdClusterSize1(i) = sqrt(sum(clusterDist1(:,i).*(1:length(EDGEScnumb)-1).^2')/sum(clusterDist1(:,i))...
                                    -Stats.meanClusterSize1(i)^2);
        Stats.stdClusterSize2(i) = sqrt(sum(clusterDist2(:,i).*(1:length(EDGEScnumb)-1).^2')/sum(clusterDist2(:,i))...
                                    -Stats.meanClusterSize2(i)^2);
        % statistics for clusters not counting individual cells (or pairs technically)                  
        Stats.meanClusterSizeGroups1(i) = sum(clusterDist1(2:end,i).*(2:length(EDGEScnumb)-1)')/sum(clusterDist1(2:end,i));
        Stats.meanClusterSizeGroups2(i) = sum(clusterDist2(2:end,i).*(2:length(EDGEScnumb)-1)')/sum(clusterDist2(2:end,i));
        Stats.stdClusterSizeGroups1(i) = sqrt(sum(clusterDist1(2:end,i).*(2:length(EDGEScnumb)-1).^2')/sum(clusterDist1(2:end,i))...
                                    -Stats.meanClusterSize1(i)^2);
        Stats.stdClusterSizeGroups2(i) = sqrt(sum(clusterDist2(2:end,i).*(2:length(EDGEScnumb)-1).^2')/sum(clusterDist2(2:end,i))...
                                    -Stats.meanClusterSize2(i)^2);

        Stats.meanClusterLength1(i) = mean([clusterLength1{:,i}]);
        Stats.meanClusterLength2(i) = mean([clusterLength2{:,i}]);
        Stats.stdClusterLength1(i) = std([clusterLength1{:,i}]);
        Stats.stdClusterLength2(i) = std([clusterLength2{:,i}]);
    end
    
    
    h3 = figure; 
    bar(1:nBars,(1:nBars).*clusterDist1(1:nBars,1)'/sum((1:params.N).*clusterDist1(:,1)'))
    title('Inital probability of cell being in cluster by size for angle 1')
    xlabel('group size')
    ylabel('probability')
    ylim([0 1])
    if SAVEFIG
        savefig(['Probability of cell being in cluster by size initial time'])
        saveas(h3,['Probability of cell being in cluster by size initial time.png'],'png')
    end
    
    h4 = figure; 
    bar(1:nBars,(1:nBars).*clusterDist1(1:nBars,end)'/sum((1:params.N).*clusterDist1(:,end)'))
    title('Final probability of cell being in cluster by size for angle 1')
    xlabel('group size')
    ylabel('probability')
    ylim([0 1])
    if SAVEFIG
        savefig(['Probability of cell being in cluster by size final time'])
        saveas(h4,['Probability of cell being in cluster by size final time.png'],'png')
    end
    
    if SAVE
        save('cluster size distribution with cluster diameters','clusterSizes1','clusterSizes2','clusterDist2','clusterDist1',...
            'Stats','clusterLength1','clusterLength2')
    end
    
    
end


function [clusterSizes1,clusterSizes2] = clusterAnalysis(identificationString,varargin)
%CLUSTERANALYSIS Take sim data and parameters and calculates clusters for
%each time step. Additionally calculates distribution of cluster sizes.


    PLOT = true;
    SAVE = true;
    TITLE = ['Clusters in final time'];
    
    % NOTE: change params.dt below to match min distance for different time steps
    stepScale = 1;
    
    Files = dir(fullfile('./',['*' identificationString '*.mat']));
    
    minN = 3; %minimum number of neightbors to count in a cluster
    
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
    clusterMajorAxis1 = cell(length(Files),length(data));
    clusterMajorAxis2 = cell(length(Files),length(data));
    EDGES = -1.5:1:(params.N/minN+0.5);
    EDGEScnumb = 0.5:1:params.N+0.5;
    
    
    for i = 1:length(data)
        
        ang = unique(data{i}.o);
        
        for k = ang'

            X = data{i}.pos(data{i}.o == k,:);
            idx = dbscan(X,params.len/1.9+params.dt*params.vbar*stepScale,minN);
            clusterSizes1{1,i} = histcounts(idx,EDGES);
            clusterSizes2{1,i} = histcounts(idx,EDGES);
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
                
                idx = dbscan(X,params.len/1.9+params.dt*params.vbar*stepScale,minN);
                if k == ang(1)
                    clusterSizes1{f,i} = histcounts(idx,EDGES);
                elseif k == ang(2)
                    clusterSizes2{f,i} = histcounts(idx,EDGES);
                end
                
                if f == length(Files) && i == length(data) && PLOT
                    figure; gscatter(X(:,1),X(:,2),idx)
                    title([TITLE ' for angle ' num2str(k)])
                end

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

%         Stats.meanClusterMajorAxis1(i) = sum(clusterMajorAxis1(2:end,i).*(2:length(EDGEScnumb)-1)')/sum(clusterDist1(2:end,i));
%         Stats.meanClusterMajorAxis2(i) = sum(clusterMajorAxis2(2:end,i).*(2:length(EDGEScnumb)-1)')/sum(clusterDist2(2:end,i));
%         Stats.stdClusterMajorAxis1(i) = sqrt(sum(clusterMajorAxis1(2:end,i).*(2:length(EDGEScnumb)-1).^2')/sum(clusterDist1(2:end,i))...
%                                     -Stats.meanClusterSize1(i)^2);
%         Stats.stdClusterMajorAxis2(i) = sqrt(sum(clusterDist2(2:end,i).*(2:length(EDGEScnumb)-1).^2')/sum(clusterDist2(2:end,i))...
%                                     -Stats.meanClusterSize2(i)^2);
    end
    
    if SAVE
        save('cluster size distribution normal dt','clusterSizes1','clusterSizes2','clusterDist2','clusterDist1',...
            'Stats')
    end
    
    
end


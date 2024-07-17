function [clusterCellCount] = clusterAnalysisShort(folder_name,dataFileName)
% Calculated number of cells in clusters and clusters themselves for given
% folder
    SAVEDATA = false;
    minN = 2; %minimum number of neightbors to count in a cluster
    
    
    load([folder_name filesep dataFileName],'dataCat');
    load([folder_name filesep 'data for run 1.mat'],'params');
    ang = unique(dataCat{1}.o);

    EDGES = -1.5:1:(params.N/minN+0.5); % edges of bins for cluster ids
    clusterCellCount = cell(length(dataCat),2);
    for j = 1:length(dataCat)
        j
        data = dataCat{j};

        runs = length(data.pos)/params.N;
        clusterSizes1 = cell(runs,1);
        clusterSizes2 = cell(runs,1);
        for i = 1:runs
            d = data.pos((i-1)*params.N+1:i*params.N,:);
            o = data.o((i-1)*params.N+1:i*params.N);
            data_slice1 = d(o == ang(1),:);
            data_slice2 = d(o == ang(2),:);

%             idx1 = dbscan(data_slice1,params.len/1.9+params.dt*params.vbar*stepScale,minN);
%             idx2 = dbscan(data_slice2,params.len/1.9+params.dt*params.vbar*stepScale,minN);
            idx1 = dbscan(data_slice1,params.len*1.05,minN);
            idx2 = dbscan(data_slice2,params.len*1.05,minN);

            clusterSizes1{i} = histcounts(idx1,EDGES);
            clusterSizes2{i} = histcounts(idx2,EDGES);


        end

        clusterCellCount{j,1} = clusterSizes1;
        clusterCellCount{j,2} = clusterSizes2;

    end
    
    if SAVEDATA
        save([folder_name filesep 'Cluster sizes'],'clusterCellCount','-v7.3')
    end


end
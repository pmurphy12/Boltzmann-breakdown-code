function [sigma_x,sigma_y] = diffusionAnalysis(folder_name,dataFileName)
% Calculated diffusion coefficient (actually the standard deviation of the
% step distance in units of length) in x- and y-directions based on
% individual paths of simulated rods. This is calculated regardless of
% orientation.

    SAVEDATA = true;
    PLOT = true;
    thresh = 300;
    
    
    load([folder_name filesep dataFileName],'dataCat');
    load([folder_name filesep 'data for run 1.mat'],'params');
    ang = unique(dataCat{1}.o);
    
    L = length(dataCat{1}.o)*length(dataCat);
    L2 = length(dataCat{1}.o);
    distDiff = zeros(L,2);
    
    for j = 1:length(dataCat)-1
        j
        data1 = dataCat{j};
        data2 = dataCat{j+1};

        distDiff(L2*(j-1)+1:L2*j,:) = data2.pos-data1.pos;

    end
    
    distDiff = distDiff(sqrt(distDiff(:,1).^2 + distDiff(:,2).^2) < thresh,:);
    
    sigma_x = std(distDiff(:,1));
    sigma_y = std(distDiff(:,2));
    counts_x = histcounts(distDiff(:,1));
    counts_y = histcounts(distDiff(:,2));
    
    if SAVEDATA
%         save([folder_name filesep 'Cluster sizes'],'clusterCellCount','-v7.3')
        save([folder_name filesep 'step size standard deviation'],'sigma_x','sigma_y','counts_x','counts_y','distDiff')
    end
    
    if PLOT
        h1 = figure;
        histogram(distDiff(:,1))
        ylabel('count')
        xlabel('x-direction step distance')
        
        h2 = figure;
        histogram(distDiff(:,2))
        ylabel('count')
        xlabel('y-direction step distance')
    end


end
function [diffSpread] = diffusionAnalysisSpread2(folder_name)
% Calculated diffusion coefficient (actually the standard deviation of the
% step distance in units of length) in x- and y-directions based on
% individual paths of simulated rods. This is calculated regardless of
% orientation.

    dataFileName = '1D data combined.mat';

    SAVEDATA = true;
    PLOT = true;
    fig_name = 'std (from profile) of left-moving profile over time';
    SAVEFIG = true;
    
    load([folder_name filesep dataFileName],'data1D');
    load([folder_name filesep 'data for run 1.mat'],'params');
    
    L1 = length(data1D{1,1});
    L2 = size(data1D,2);
    diffSpread = zeros(L2,1);
    x = linspace(0,params.lengthGx,L1);
    time = linspace(0,params.tend,L2);
    
    for j = 1:L2
        j
        data = data1D{1,j};
%         data = data - min(data(2:end-1)); % subtract background density
%         data(1) = 0;
%         data(end) = 0;
%         if j == 18
%             5+1;
%         end
%         [~,~,diffSpread(j),~] = circStats(x,data,params.lengthGx);
        [~,~,diffSpread(j)] = eucStats(x,data);
        

    end
    
    
    if SAVEDATA
%         save([folder_name filesep 'Cluster sizes'],'clusterCellCount','-v7.3')
        save([folder_name filesep 'step size standard deviation from profile'],'diffSpread','time','x')
    end
    
    if PLOT
        h = figure;
        plot(time,diffSpread)
        xlabel('time')
        ylabel('Standard Deviation')
    end
    
    if SAVEFIG
        savefig(h,[folder_name filesep fig_name])
        saveas(h,[folder_name filesep fig_name '.png'],'png')
        exportgraphics(h,[folder_name filesep fig_name '.eps'],'ContentType','vector')
    end


end
function [diffSpread] = diffusionAnalysisSpread(folder_name,dataFileName)
% Calculated diffusion coefficient (actually the standard deviation of the
% step distance in units of length) in x- and y-directions based on
% individual paths of simulated rods. This is calculated regardless of
% orientation.

    SAVEDATA = true;
    PLOT = true;
    fig_name = 'std of left-moving profile over time';
    SAVEFIG = true;
    
    load([folder_name filesep dataFileName],'dataCat');
    load([folder_name filesep 'data for run 1.mat'],'params');
    ang = unique(dataCat{1}.o);
    
    L = length(dataCat);
    diffSpread = zeros(L,1);
    time = zeros(L,1);
    
    for j = 1:L
        j
        data = dataCat{j};
        pos = data.pos(data.o == ang(1),1);

        diffSpread(j) = std(pos);
        time(j) = data.t;

    end
    
    
    if SAVEDATA
%         save([folder_name filesep 'Cluster sizes'],'clusterCellCount','-v7.3')
        save([folder_name filesep 'step size standard deviation'],'diffSpread','time')
    end
    
    if PLOT
        h = figure;
        plot(time,diffSpread)
        xlabel('time')
        ylabel('left-moving profile std')
    end
    
    if SAVEFIG
        savefig(h,[folder_name filesep fig_name])
        saveas(h,[folder_name filesep fig_name '.png'],'png')
        exportgraphics(h,[folder_name filesep fig_name '.eps'],'ContentType','vector')
    end


end
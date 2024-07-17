function [] = plotWave(folder_name)
% generates movie of left-and right-moving profiles

%% setup
    dataFileName = '1D data combined.mat';
    
    load([folder_name filesep dataFileName],'data1D');
    load([folder_name filesep 'data for run 1.mat'],'params');
    
    L1 = length(data1D{1,1});
    Lu = max([data1D{:}]);
    x = linspace(0,params.lengthGx,L1);
    
    curr_folder = pwd;
    
    %% 1D movies
    
    cd(folder_name)
    
    v = VideoWriter(['wave movie']);
    v.FrameRate = 5;
    open(v)
    
    figure
    for i=1:size(data1D,2)
        plot(x,data1D{1,i})
        hold on
        plot(x,data1D{2,i})
        hold off
        axis([0 x(end) 0 1.1*Lu])
        xlabel('x')
        ylabel('Density')
        M = getframe(gcf);
        writeVideo(v,M)
    end
    close(v)

    cd('../')


end


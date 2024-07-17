function [cfh] = angleHist(tracksData,NBins,movie,save_folder,save_name)
% angleHist creates a movie of the histogram of cell orientations averaged
% over space.
%   tracksData: a cell array where each cell is a structure containing cell
%   position (pos) and orientation data (o) at a certain time (t).
    
    
    if nargin < 3
        movie = false;
    end
    
    if ~isempty(save_folder)
        curr_folder = pwd;
        cd(save_folder)
    end
    
    if movie    
        v = VideoWriter(strcat('orientation histogram ',save_name));
        v.FrameRate = 10;
        open(v)
    end
    
    figure;
    for i = 1:length(tracksData)
       chf = histogram(tracksData{i}.o,NBins,'Normalization','probability');
       axis([-pi-0.2 pi+0.2 0 0.1]);
       xlabel('Orientation');
       ylabel('Cell Count');
       title('Orientation Histogram')
       drawnow
       if movie
            M = getframe(gcf);
            writeVideo(v,M)
       end
    end
    
    if movie
        close(v)  
    end
                
    if ~isempty(save_folder)
        cd(curr_folder)
    end 
    
    
end


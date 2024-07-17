function [ratioData] = coneRatioDiff(tracksData,params,r,movie,save_folder,save_name)
% coneRatio checks the ratio of densities within a radius r of location (x,y) 
% between antipodal orientations is independent (or asymptotically
% independent) of angle (i.e. rho+/rho- is a function of only x, y, and t).

%   tracksData: a cell array where each cell is a structure containing cell
%   position (pos) and orientation data (o) at a certain time (t).
%
%   params: a structure containing parameters used in paralellogram
%   collision model, particularly grid size and grid spacing.
%
%   r: spatial radius over which to compute ratio (0.1?).
    
    
    if nargin < 3
        movie = false;
    end
    
    if ~isempty(save_folder)
        curr_folder = pwd;
        cd(save_folder)
    end
    
    if movie    
        v = VideoWriter(strcat('ratio histogram ',save_name));
        v.FrameRate = 10;
        open(v)
    end
    
    NBins = 3;
    NGrid = 10;
    angBins = linspace(0,pi,NBins+1);
    ratioData = cell(1,length(tracksData));
    [X,Y] = meshgrid(linspace(0,params.grid*params.lengthG,NGrid),linspace(0,params.grid*params.lengthG,NGrid));
    
    figure('Position',[100 50 900 720]);
    for i = 1:length(tracksData)
        ratio = zeros(NGrid,NGrid,NBins);
        searchT = createns(tracksData{i}.pos);
        for j = 1:numel(X)
            idx = rangesearch(searchT,[X(j),Y(j)],r,'SortIndices',false);
%             ori = mod(tracksData{i}.o(idx{1})+pi-0.6,2*pi)-pi;
            ori = tracksData{i}.o(idx{1});
            
            for k = 1:length(angBins)-1
                [I,J] = ind2sub(size(X),j);
%                 ratio(I,J,k) = sum(ori >= angBins(k) & ori < angBins(k+1))...
%                                 /sum(ori >= angBins(k)-pi & ori < angBins(k+1)-pi);
                ratio(I,J,k) = sum(ori >= angBins(k) & ori < angBins(k+1))...
                                /(sum(ori >= angBins(k)-pi & ori < angBins(k+1)-pi)...
                                +sum(ori >= angBins(k) & ori < angBins(k+1)));
            end
        end
        
        ratioData{i} = ratio;
        
        for ii = 1:NBins
            subplot(2,2,ii);
            bar3(ratio(:,:,ii))
            axis([-1 NGrid+1 -1 NGrid+1 0 1]);
            xlabel('x');
            ylabel('y');
            zlabel(['Ratio for orientation bin ' num2str(angBins(ii)) '-' num2str(angBins(ii+1))])
            title('Ratio as a function of angle')
        end
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


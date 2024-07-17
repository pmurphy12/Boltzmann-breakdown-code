function  MC_movieWriter(PosData,OrientData,NColorBins,MovieNumber,MovieTitle,grid,lengthGxy,len,varargin)
% Generates movies for alignmentMonteCarlo scripts by binning particles into
% NColorBins bins based on orientation, then plotting each individually and
% saving the frame to a movie. It as plots all colors at once and saves
% frame into a movie of all particles. 
%
% Orientation and position data are fed in through OrientData and PosData
% respectively. Both are cell arrays size Nt x 1, where Nt is the number of
% time steps. Each cell is N x 2, where N is the number of particles.
%
% MovieNumber indexes movies accross individual runs.
    
    currDir = pwd;
    if length(varargin) == 2 && exist([varargin{2}],'dir') == 7
        changeDir = varargin{2};
        disp(changeDir)
        cd(changeDir)
    elseif exist([pwd '/Results'],'dir') == 7
        cd './Results'
    end
    
    dims = size(PosData);
    Nt = max(dims,[],'all');
    
    if ~isempty(varargin)
        NColorBins = length(varargin{1})-1;
    end
    
   % Create movie names and bin data
    for j = 1:NColorBins
        moviename{j} = [MovieTitle sprintf('movie AngBins_%i_%d', MovieNumber, j)];
        col{j} = color_chooser(j);
    end
    
    PosBinnedAll = cell(Nt,1); %Preallocate
    OrientBinnedAll = cell(Nt,1);
    
    % Check for optional arguments
    if ~isempty(varargin)
        for i = 1:Nt
            [PosBinnedAll{i},OrientBinnedAll{i}] = GenColorBins(PosData{i},...
                OrientData{i},NColorBins,varargin{1});
        end
    elseif isempty(varargin)
        for i = 1:Nt
            [PosBinnedAll{i},OrientBinnedAll{i}] = GenColorBins(PosData{i},...
                OrientData{i},NColorBins);
        end
    else
        error('Exceeded the maximum number of inputs.')
    end
    
    % Create individual movies for each bin.
    for j = 1:NColorBins
        
        v = VideoWriter(moviename{j});
        open(v)
        
        figure('Renderer', 'painters', 'Position', [10 10 900 900*lengthGxy(2)/lengthGxy(1)])
        
        p = Progress(Nt)
        for i = 1:Nt
            p.d(i) % Progress bar update
            
            % Plot vectors representing each cell in given bin.
            quiver(PosBinnedAll{i}{j}(:,1)-len/2*OrientBinnedAll{i}{j}(:,1),...
                PosBinnedAll{i}{j}(:,2)-len/2*OrientBinnedAll{i}{j}(:,2),...
                len*OrientBinnedAll{i}{j}(:,1),len*OrientBinnedAll{i}{j}(:,2),...
                0,'Color',col{j})
            
            % Set axis
            ax = gca;
            ax.YLim = [0,grid*lengthGxy(2)];
            ax.XLim = [0,grid*lengthGxy(1)];
            drawnow         %Updates figure within for loop.
            
            % Get frame for movie
            M = getframe(gcf);
            writeVideo(v,M)
        end
        % Close movie object
        close(v)
        p.done()
    end
    
    % Create full movie
    v = VideoWriter([MovieTitle sprintf('movie full %d', MovieNumber)]);
    open(v)
    
    figure('Renderer', 'painters', 'Position', [10 10 900 900*lengthGxy(2)/lengthGxy(1)])
    
    p = Progress(Nt)
    for i=1:Nt
        p.d(i) % Progress bar update
        
        % Plot vectors representing each cell in given bin.
        for j = 1:NColorBins
            quiver(PosBinnedAll{i}{j}(:,1)-len/2*OrientBinnedAll{i}{j}(:,1),...
                PosBinnedAll{i}{j}(:,2)-len/2*OrientBinnedAll{i}{j}(:,2),...
                len*OrientBinnedAll{i}{j}(:,1),len*OrientBinnedAll{i}{j}(:,2),...
                0,'Color',col{j})
            hold on
        end
        hold off
        
        % Set axis
        ax = gca;
        ax.YLim = [0,grid*lengthGxy(2)];
        ax.XLim = [0,grid*lengthGxy(1)];
        drawnow         %Updates figure within for loop.
        
        % Get frame for movie
        M = getframe(gcf);
        writeVideo(v,M)
    end
    % Close movie object
    close(v)
    p.done()
    
     cd(currDir);
end


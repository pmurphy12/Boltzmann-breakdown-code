% Code for Monte Carlo simulations of cell orientiation alignment through
% interactions. Assumes periodic domain. Calculates and stores data such as histogram of
% orientations as function of time as well as spatial distribution of
% cells. Records movie if option selected.

% Aligment is assumed to occur instantly. If cell i comes withing a certain
% distance of cell j, the orientation angle of cell i is set to that of
% cell j. Cell midpoints are assumed to be at coordinates (x_i,y_i) with
% orientation vector \xi_i.
% Based on a more detailed interaction using parallelogram vector calculation.
% coll_type is 'Polar' or 'Nematic'.
%
% Uses exact forward-looking collision scheme that tests for
% collisions before updating position,reorienting, then updating position again.
%
% Each agent cell has an individual refractory period given by p_refract.
% Cell reorientation is assumed to occur at the head of the cell.

%SORTINDICIES DURING KNN SEARCH CURRENTLY SET TO FALSE


function [params,tracksData,collisions] = alignMC_quasirefractory(MovieNumber,len,N,L,tend,m,IC_type,coll_type,varargin)

    p = inputParser;
    addOptional(p,'InitialMean',[],@isnumeric);
    addOptional(p,'InitialSD',[],@isnumeric)
    addOptional(p,'InitialProp',[],@isnumeric);
    addOptional(p,'SaveFolder','',@ischar);
    addOptional(p,'InitialData',{},@isstruct);
    addOptional(p,'BinLoc',[],@isnumeric);
    addOptional(p,'RunNumb',[],@isnumeric);
    addOptional(p,'MovieTitle','',@ischar);
    addOptional(p,'yShift',false,@islogical);
    addOptional(p,'yShift_supp',0,@isnumeric);
    addOptional(p,'yShift_all',false,@islogical);
    addOptional(p,'NoiseLevel',0,@isnumeric);
    addOptional(p,'NoCollisions',false,@islogical);
    
    parse(p,varargin{:});  
    
    rng('shuffle')
%     rng('default')
    
    parallel.defaultClusterProfile
    parallel.defaultClusterProfile('local')
    %c = parcluster
    %poolobj = parpool(c)
        
    Plotting = 0; % Turn on/off regular plotting.
    ColorBinsPlot = 1; % Turn on/off color bin plotting. Need Plotting=1 as well.
    
    CMovie = 0; % Turn on/off movie making for color bins.
    NColorBins = 2; % Number of orientation bins for plotting with color.
    
    Collect = 1; % Collect data for mean, variance, as well as histogram data.
    

    
    params.EPS = 1e-8; % Tolerance for rounding errors.
    params.grid = 1; % Grid number.
    params.lengthGx = 400; % Grid length x-direction, dimensional.
    params.lengthGy = L; % Grid length y-direction, dimensionlal.
    
    params.N = N; %14400; % Number of cells in simulation.
    params.len = len; % Cell length in micrometers. Sets interaction distance.
    params.vbar = sqrt(2); %10, Cellular velocity. Micrometers per minute.
    
    params.p_refract = params.len/params.vbar;
    
    params.tend = tend; % Simulation length in minutes.
    params.m = m;
    params.dt = params.len/params.vbar/(10*params.m);%params.len/params.vbar/2;%0.0025; % time step size
    params.Nt = ceil(params.tend/params.dt); 
    
    params.K = params.N*params.len^2/params.lengthGy^2/2 % Interaction coefficient kappa in kinetic model for reference.
    
    params.sigma = p.Results.NoiseLevel; % Standard deviation for gaussian noise in orientation
    params.yshift = p.Results.yShift_supp;
    
    % Preallocate for data collection
    Sample = min(100,floor(params.Nt/1)); %params.Nt;% How many samples to take over the simulation time.
    SampleRate = params.tend/Sample;
    SampleNum = 1;
    MSample = min(300,floor(params.Nt/1)); % How many frame to collect over the simulation time.
    MSampleRate = params.tend/MSample;
    MSampleNum = 1;
    
    params.KNNsteps = floor(40*m); %40
    
    params.saveInterval = floor(100*m); %100
    if isempty(p.Results.RunNumb)
        RunNumb = 1;
    else
        RunNumb = p.Results.RunNumb+1;
    end
        
    if ~isempty(p.Results.SaveFolder)
        save_state = true;
    else
	save_state = false;
    end

    % Set data bins
%     Nangle = 36; % Number of bins for angle data. 36 gives 10 degree bin size.
%     BinPos = linspace(-pi,pi,Nangle); % Bin center locations.
    
    if Collect
        tracksData = cell(1,Sample+1);
    end
    
    % Positional data binned for crude color binning. Feeds into
    % MC_movieWriter
    if CMovie
        PosData = cell(MSample+1,1);
        OrientData = cell(MSample+1,1);
    end


    if isempty(p.Results.InitialData)
        % Initialize positions. 'Uniform' 'Test'
        [Pos,Ang,BinLocations] = initializeCells(params,IC_type,...
            'InitialMean',p.Results.InitialMean,'InitialSD',p.Results.InitialSD,...
            'InitialProp',p.Results.InitialProp);
        Orient = [cos(Ang) sin(Ang)];
%         collisions = zeros(params.N,1);
        collisions = sparse(params.N,params.Nt);
        refract_times = params.p_refract*ones(params.N,1);
    else
        Pos = p.Results.InitialData.pos;
        Orient = [cos(p.Results.InitialData.o) sin(p.Results.InitialData.o)];
        BinLocations = p.Results.BinLoc;
        Ang = atan2(Orient(:,2),Orient(:,1));
        collisions = p.Results.InitialData.collisions;
        refract_times = p.Results.InitialData.refract_times;
    end
    
    if Plotting
        plotfig = figure;
        hAxes = gca;
        [PosBinned,OrientBinned] = GenColorBins(Pos,Orient,NColorBins,BinLocations);
        
        for j = 1:length(PosBinned)
            
            quiver(hAxes,PosBinned{j}(:,1)-params.len*OrientBinned{j}(:,1),PosBinned{j}(:,2)-params.len*OrientBinned{j}(:,2),...
                params.len*OrientBinned{j}(:,1),params.len*OrientBinned{j}(:,2),0,'Color',color_chooser(j))
            hold(hAxes,'on')
            
        end
        hold(hAxes,'off')

        hAxes.YLim = [0,params.grid*params.lengthGy];
        hAxes.XLim = [0,params.grid*params.lengthGx];
%         axis([0 params.lengthGx 0 params.lengthGy])
        drawnow 
        pause
    end
    
    if Collect
        % Collect initial orientation and position data.
        tracksData{SampleNum}.pos = Pos;
        tracksData{SampleNum}.o = Ang;
        tracksData{SampleNum}.t = 0;
    end
    
    if CMovie
        % Calculate initial orientation data.
        PosData{MSampleNum} = Pos;
        OrientData{MSampleNum} = Orient;
    end
    

    pr = Progress(params.Nt) % Starts progress bar in Command Window.
    
    % Main for loop
    
    for i = 1:params.Nt
        pr.d(i) % Progress bar update
        
        tCurr = i*params.dt; % Simulation time tracker.
        
        I_refract = refract_times < params.p_refract;
%         tempPos = Pos; % FOR DEBUG
        
        % Alternate collision check using createns and rangesearch
        % (non-exhaustive). Catches collisions in the middle of movement.
        if ~p.Results.NoCollisions
            if mod(i,params.KNNsteps) == 1
                r = 2*params.len+2*params.KNNsteps*params.vbar*params.dt;
                Lx = params.lengthGx*params.grid;
                Ly = params.lengthGy*params.grid;
                % Check which cells are near boundaries
                yI = Pos(:,2) < r | Pos(:,2) > params.grid*params.lengthGy-r;
                xI = Pos(:,1) < r | Pos(:,1) > params.grid*params.lengthGx-r;
                posMdl_1 = createns(Pos,'NSMethod','kdtree');
                posMdl_2 = createns([mod(Pos(xI,1)+Lx/2,Lx),Pos(xI,2)],'NSMethod','kdtree');
                posMdl_3 = createns([Pos(yI,1),mod(Pos(yI,2)+Ly/2,Ly)],'NSMethod','kdtree');
                posMdl_4 = createns([mod(Pos(xI & yI,1)+Lx/2,Lx),mod(Pos(xI & yI,2)+Ly/2,Ly)],'NSMethod','kdtree');


                [idx1,D1] = rangesearch(posMdl_1,Pos,r,'SortIndices',false);
                [idx2,D2] = rangesearch(posMdl_2,[mod(Pos(xI,1)+Lx/2,Lx),Pos(xI,2)],r,'SortIndices',false);
                [idx3,D3] = rangesearch(posMdl_3,[Pos(yI,1),mod(Pos(yI,2)+Ly/2,Ly)],r,'SortIndices',false);
                [idx4,D4] = rangesearch(posMdl_4,[mod(Pos(xI & yI,1)+Lx/2,Lx),mod(Pos(xI & yI,2)+Ly/2,Ly)],r,'SortIndices',false);
%                 [idx1,D1] = rangesearch(posMdl_1,Pos,r);
%                 [idx2,D2] = rangesearch(posMdl_2,[mod(Pos(xI,1)+Lx/2,Lx),Pos(xI,2)],r);
%                 [idx3,D3] = rangesearch(posMdl_3,[Pos(yI,1),mod(Pos(yI,2)+Ly/2,Ly)],r);
%                 [idx4,D4] = rangesearch(posMdl_4,[mod(Pos(xI & yI,1)+Lx/2,Lx),mod(Pos(xI & yI,2)+Ly/2,Ly)],r);
                
                idxF = cell(length(idx1),1);

                % Create cell arrays with right lengths by filling in with empty
                % cells
                tempIdx = cell(length(idx1),1);
                tempIdx(xI) = idx2;
                idx2 = tempIdx;
                tempD = cell(length(D1),1);
                tempD(xI) = D2;
                D2 = tempD;

                tempIdx = cell(length(idx1),1);
                tempIdx(yI) = idx3;
                idx3 = tempIdx;
                tempD = cell(length(D1),1);
                tempD(yI) = D3;
                D3 = tempD;

                tempIdx = cell(length(idx1),1);
                tempIdx(xI & yI) = idx4;
                idx4 = tempIdx;
                tempD = cell(length(D1),1);
                tempD(xI & yI) = D4;
                D4 = tempD;
            end

    %         if i == 100
    %             5+1
    %         end
            
            % remove refractory cells from colliding with another cell
            idx1(I_refract) = cell(sum(I_refract),1);
            
            Ind = cell(length(idx1),1);
            Beta = cell(length(idx1),1);

            parfor ii = 1:length(idx1)

                [idx] = union(idx1{ii},idx2{ii});
                [idx] = union(idx,idx3{ii});
                [idx] = union(idx,idx4{ii});

                I = idx ~= ii;
                idx = idx(I);

                O = Orient(ii,:);
                Otest = Orient(idx,:) - O ~= 0;
                idx = idx(logical(sum(Otest,2)));


                if ~isempty(idx)

                    %  Get orientations of cells possibly colliding
                    SetOfOrient = Orient(idx,:);

                    % Get relevant spatial data of cells possibly colliding
                    P = Pos(ii,:);
                    SetOfPos = Pos(idx,:);

                    DiffInPos = min(abs(P-SetOfPos),...
                    params.grid*[params.lengthGx,params.lengthGy]-abs(P-SetOfPos));
                    T1 = P < SetOfPos;
                    T2 = abs(P-SetOfPos) > params.grid*[params.lengthGx,params.lengthGy]-abs(P-SetOfPos);
                    DiffInPos = DiffInPos.*(1-2*T1).*(1-2*T2);

                    ProjBasis = dot(repmat([-O(2),O(1)],[length(idx),1]),SetOfOrient,2);
                    Proj1 = dot(DiffInPos,repmat([-O(2),O(1)],[length(idx),1]),2);
                    Proj2 = dot(DiffInPos,[-SetOfOrient(:,2),SetOfOrient(:,1)],2);
                    
                    ParamAlpha = Proj1./ProjBasis;
                    ParamBeta = Proj2./ProjBasis;
                    
                    ParamTestA = (ParamAlpha < params.vbar*params.dt).*(ParamAlpha > -params.len); % Test for cell i hitting cell j
                    ParamTestB = (ParamBeta < params.vbar*params.dt).*(ParamBeta > 0);
                    ParamTestC = ParamAlpha < ParamBeta; % Test for if cell i will hit cell j instead of the other way around. Combined with previous tests A and B, this condition is enough to cover all cases.

                    Test = logical(ParamTestA.*ParamTestB.*ParamTestC);

                    idx = idx(Test);
                    Beta{ii} = ParamBeta(Test);
                    if length(idx) > 1
                        [Beta{ii},Imin] = min(Beta{ii});
                        idx = idx(Imin);
                    end

                    idxF{ii} = idx;
                    Ind{ii} = ii*ones(1,length(idx));

                else
                    idxF{ii} = idx';
                    Beta{ii} = [];
                    Ind{ii} = ii*ones(1,length(idx));
                end


            end

            % Get rid of collisions that cause pairs of cells to switch
            % orientations.
            I = cat(2,[Ind{:}]);
            
        else
            I = [];
        end
        
            if ~isempty(I)
                J = cat(2,[idxF{:}]);
                B = cat(2,[Beta{:}]);

                % Update CellNum1 and CellNum2 based on results of Test1 and
                % multicollision check, using updated matrix I.
                CellNum1 = I;
                CellNum2 = J;


                % Update cell positions.
                notI = setdiff(1:params.N,I);
                Pos(notI,:) = Pos(notI,:)+params.vbar*params.dt*Orient(notI,:);
                Pos(I,:) = Pos(I,:)+(params.vbar*params.dt/2).*Orient(I,:)+(params.vbar*params.dt/2).*Orient(J,:); % average collision update step
                % Pos(I,:) = Pos(I,:)+(B').*Orient(I,:)+(params.vbar*params.dt-(B')).*Orient(J,:); % Exact collisions update step

                % shift in y-direction of yShift set to true. Tests
                % clustering effects.
                if p.Results.yShift
                    Pos(I,2) = Pos(I,2) + p.Results.yShift_supp*params.len*(rand(length(Pos(I,2)),1)-1.5/5); % -1.5/5
                elseif p.Results.yShift_all
%                     Pos(:,2) = Pos(:,2) + p.Results.yShift_supp*params.len*(rand(length(Pos(:,2)),1)-1.5/5);
                    Pos(:,2) = Pos(:,2) + p.Results.yShift_supp*(randn(length(Pos(:,2)),1)); % -1.5/5
                end

                % Implement periodic boundary conditions
                Pos = Pos + (Pos<0).*params.grid.*[params.lengthGx,params.lengthGy] ...
                    - (Pos>params.grid.*[params.lengthGx,params.lengthGy]).*params.grid.*[params.lengthGx,params.lengthGy];
                
%                 diffPos = Pos(:,2)-tempPos(:,2); % FOR DEBUG
%                 diffPos = diffPos(sqrt(diffPos.^2 + diffPos.^2) < 300); % FOR DEBUG
%                 if sum(abs(diffPos) <= 0.1) % FOR DEBUG
%                     5+1; % FOR DEBUG
%                 end % FOR DEBUG

                % Update collision counts
%                 collisions(CellNum1) = collisions(CellNum1)+1;
                collisions(CellNum1,i) = 1;

                % Change orientations for cells that collided.

                %Nematic
                if strcmp(coll_type,'Nematic')
                    Ang = atan2(Orient(:,2),Orient(:,1));
                    angDiff = min(abs(Ang(CellNum1)-Ang(CellNum2)),...
                        2*pi-abs(Ang(CellNum1)-Ang(CellNum2)));
                    % Collision of cell 1 into cell 2.
                    OrientNum1 = Orient(CellNum1,:);
                    OrientNum2 = Orient(CellNum2,:);
                    OrientNum1((angDiff < pi/2),:) = OrientNum2((angDiff < pi/2),:);
                    OrientNum1((angDiff > pi/2),:) = [-OrientNum2((angDiff > pi/2),1), ...
                                                                        -OrientNum2((angDiff > pi/2),2)];
        %             if sum(angDiff <= pi/2+params.EPS & angDiff >= pi/2-params.EPS)
                    if sum(angDiff <= pi/2 & angDiff >= pi/2)
                        b = logical(binornd(1,0.5,size(OrientNum1,1),1));
                        OrientNum1((angDiff == pi/2) & b,:) = OrientNum2((angDiff == pi/2) & b,:);
                        OrientNum1((angDiff == pi/2) & ~b,:) = [-OrientNum2((angDiff == pi/2) & ~b,2), ...
                                                                        -OrientNum2((angDiff == pi/2) & ~b,1)];
                    end
                    Orient(CellNum1,:) = OrientNum1;

                % Polar
                elseif strcmp(coll_type,'Polar')
                    Orient(CellNum1,:) = Orient(CellNum2,:); % Collision of
                    % cell 1 into cell 2.

                end

                % Update refractory periods
                refract_times(CellNum1) = 0;
                refract_times = refract_times+params.dt;

                %update orientations with noise
                eps = params.sigma*randn(params.N,1);
                Ang = atan2(Orient(:,2),Orient(:,1)) + eps;
                Orient = [cos(Ang), sin(Ang)];
                

            else
                % Update Positions
                Pos = Pos+params.vbar*params.dt*Orient;
                
                if p.Results.yShift_all
%                     Pos(:,2) = Pos(:,2) + p.Results.yShift_supp*params.len*(rand(length(Pos(:,2)),1)-1.5/5);
                    Pos(:,2) = Pos(:,2) + p.Results.yShift_supp*(randn(length(Pos(:,2)),1)); % -1.5/5
                end
                
                % Implement periodic boundary conditions
                Pos = Pos + (Pos<0).*params.grid.*[params.lengthGx,params.lengthGy] ...
                    - (Pos>params.grid.*[params.lengthGx,params.lengthGy]).*params.grid.*[params.lengthGx,params.lengthGy];
                
                % Update refractory periods
                refract_times = refract_times+params.dt;

                %update orientations with noise
                eps = params.sigma*randn(params.N,1);
                Ang = atan2(Orient(:,2),Orient(:,1)) + eps;
                Orient = [cos(Ang), sin(Ang)];

            end
            
        
        

        
        % Collect data
        if tCurr >= SampleNum*SampleRate && Collect
            
            Ang = atan2(Orient(:,2),Orient(:,1));
            
            SampleNum = SampleNum+1; % Update sample number
            tracksData{SampleNum}.pos = Pos;
            tracksData{SampleNum}.o = Ang;
            tracksData{SampleNum}.t = tCurr;

        end

        if CMovie && (tCurr >= MSampleNum*MSampleRate)
            MSampleNum = MSampleNum+1; % Update sample number
            PosData{MSampleNum} = Pos;
            OrientData{MSampleNum} = Orient;
        end
        
        % Save data if save fodler specified
        if save_state && mod(i,params.saveInterval) == 0
            save([p.Results.SaveFolder filesep 'savestate' num2str(RunNumb)],'params','i','tracksData','SampleNum','BinLocations','RunNumb','collisions')
        end
        
        

        % Regular Plotting
        if Plotting && ~ColorBinsPlot % && tCurr >= SampleNum*SampleRate
            % Vector plot using cell midpoint and orientation vector scaled
            % to params.len/2.
%             quiver(Pos(:,1)-params.len/2*Orient(:,1),Pos(:,2)-params.len/2*Orient(:,2),...
%                 params.len*Orient(:,1),params.len*Orient(:,2),0)

            % Plotting two populations in different colors.
            [PosBinned,OrientBinned] = GenColorBins(Pos,Orient,1);
            
            quiver(hAxes,PosBinned{1}(:,1)-params.len*OrientBinned{1}(:,1),PosBinned{1}(:,2)-params.len*OrientBinned{1}(:,2),...
                params.len*OrientBinned{1}(:,1),params.len*OrientBinned{1}(:,2),0,'Color',[0.8500 0.3250 0.0980])
            
            hAxes.YLim = [0,params.grid*params.lengthGy];
            hAxes.XLim = [0,params.grid*params.lengthGx];
            drawnow
        end
        
        
            % Color Plotting
        if Plotting && ColorBinsPlot % && tCurr >= SampleNum*SampleRate
            % Plot each orientation bin with a different color in a different figure.
            
            
            % Orientation binning based on NColorBins. Bin orientations into 
            % NColorBins bins.
            [PosBinned,OrientBinned] = GenColorBins(Pos,Orient,NColorBins,BinLocations);
            
            for j = 1:length(PosBinned)

                % Vector plot using cell midpoint and orientation vector scaled
                % to params.len/2.
                quiver(hAxes,PosBinned{j}(:,1)-params.len*OrientBinned{j}(:,1),PosBinned{j}(:,2)-params.len*OrientBinned{j}(:,2),...
                    params.len*OrientBinned{j}(:,1),params.len*OrientBinned{j}(:,2),0,'Color',color_chooser(j))
                hold(hAxes,'on')
            end      
            hold(hAxes,'off')
            hAxes.YLim = [0,params.grid*params.lengthGy];
            hAxes.XLim = [0,params.grid*params.lengthGx];
            drawnow
            pause(0.5)
        end
        
        

        
        
        
    end
     pr.done() % Carrage return for progress bar display in Command Window.
    
    
    % Saving data
%     save('fileName','variableName')
    
    
    if CMovie && exist('BinLocations','var')
        if save_state
            MC_movieWriter(PosData,OrientData,NColorBins,MovieNumber,p.Results.MovieTitle,params.grid,[params.lengthGx,params.lengthGy],params.len,BinLocations,p.Results.SaveFolder)
        else
            MC_movieWriter(PosData,OrientData,NColorBins,MovieNumber,p.Results.MovieTitle,params.grid,[params.lengthGx,params.lengthGy],params.len,BinLocations)
        end
    elseif CMovie
        MC_movieWriter(PosData,OrientData,NColorBins,MovieNumber,p.Results.MovieTitle,params.grid,[params.lengthGx,params.lengthGy],params.len)
    end
    
    


    
end




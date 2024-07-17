function [data2D,data1D,X,Y] = MCkernelPlot_2d_to_1d(data,params,varargin)
%MCkernelPlot uses data from full Monte Carlo
%simulations (see alignmentMonteCarloParallelogram.m). It approxiamtes the
%probability distribution f(x,y,\xi,t) using a Gaussian kernel for position
%and orientation data for each cell simulated in the  Monte Carlo
%simulation.
%   data: a 1 x tsteps cell array containing x and y position data, 
%   orientation data, time and id [data for each cell in simulation.
%   params: a structure containing parameters used in the MC simulation.
    
    p = inputParser;
    
    
    %addOptional(p,'DistanceCutoff', 40, @isnumeric);
%     addOptional(p,'Debug', 0, @islogical);
%     addOptional(p,'StandardDevX', sqrt(0.007), @isnumeric); % for using exp func %sqrt(0.005)
%     addOptional(p,'StandardDevX', 0.005, @isnumeric); % for using mvnpdf
%     addOptional(p,'StandardDevO', 1, @isnumeric);
    addOptional(p,'LengthX', 1, @isnumeric);
    addOptional(p,'LengthY',1,@isnumeric);
    addOptional(p,'save_name','',@ischar);
    addOptional(p,'simType',' Trianglewave Test',@ischar);
    addRequired(p,'MCdata');
    addRequired(p,'MCparams');
    
    parse(p,data,params, varargin{:} );
    
    save_name = p.Results.save_name;
    simType = p.Results.simType;
    DATE = char(datetime('today','Format','MM-dd-yyyy'));
    
    N = 2^8;
    Lx = params.grid*params.lengthGx;
    Ly = params.grid*params.lengthGy;
%     dt = params.dt;
    bandwidth = 0.01;
    bwd = ([bandwidth bandwidth] ./ [Lx Ly]).^2;
    
    data2D = cell(2,length(data));
    data1D = cell(2,length(data));
    
    pr = Progress(length(data)) % Starts progress bar in Command Window.
    for i = 1:length(data)
        pr.d(i) % Progress bar update
%         i
        Pos = data{i}.pos;
        Ang = data{i}.o;
        ang = unique(Ang);
%         tCurr = data{i}.t;
        pos1 = Pos(Ang == ang(1),:);
        pos2 = Pos(Ang == ang(2),:);
        
        [~,data2D{1,i},X,Y] = kde2d(pos1,N,[0 0],[Lx Ly],bwd(1),bwd(2));
        [~,data2D{2,i},~,~] = kde2d(pos2,N,[0 0],[Lx Ly],bwd(1),bwd(2));
        data2D{1,i} = data2D{1,i}*N^2/2;
        data2D{2,i} = data2D{2,i}*N^2/2;
        
%         PDF_1d{1,i} = trapz(Y(:,1),PDF{1,i},1)/Ly;
%         PDF_1d{2,i} = trapz(Y(:,1),PDF{2,i},1)/Ly;
        
        data1D{1,i} = ksdensity(pos1(:,1),X(1,:))/2;
        data1D{2,i} = ksdensity(pos2(:,1),X(1,:))/2;
%         [data1D{1,i},~,bw1] = ksdensity(pos1(:,1),X(1,:),'Bandwidth',0.007);
%         data1D{1,i} = data1D{1,i}/2;
%         data1D{2,i} = ksdensity(pos2(:,1),X(1,:),'Bandwidth',0.007)/2;
%         bw1

    end
    pr.done() % Carrage return for progress bar display in Command Window.
    
    save([save_name simType ' after KDE ' DATE],'data2D','data1D','X','Y')
    
    %% 2D movies
%     v = VideoWriter([save_name simType ' 2D after kde2d ' DATE ' rho_1 video 75pi-2']);
%     v.FrameRate = 5;
%     open(v)
%     
%     figure
%     for i=1:size(data2D,2)
%         mesh(X,Y,data2D{1,i})
%         axis([0 Lx 0 Ly 0 max(data2D{1,1},[],'all')])
%         xlabel('x')
%         ylabel('y')
%         zlabel('Density')
%         M = getframe(gcf);
%         writeVideo(v,M)
%     end
%     close(v)
%     
%     v = VideoWriter([save_name simType ' 2D after kde2d ' DATE ' rho_2 video 25pi-2']);
%     v.FrameRate = 5;
%     open(v)
%     
%     figure
%     for i=1:size(data2D,2)
%         mesh(X,Y,data2D{2,i})
%         axis([0 Lx 0 Ly 0 max(data2D{2,1},[],'all')])
%         xlabel('x')
%         ylabel('y')
%         zlabel('Density')
%         M = getframe(gcf);
%         writeVideo(v,M)
%     end
%     close(v)
    
    %% 1D movies
    v = VideoWriter([save_name simType ' KDE ' DATE ' rho_1 1D']);
    v.FrameRate = 5;
    open(v)
    
    figure
    for i=1:size(data1D,2)
        plot(X(1,:),data1D{1,i})
        axis([0 Lx 0 max(union(data1D{1,1},data1D{1,end}),[],'all')])
        xlabel('x')
        ylabel('Density')
        M = getframe(gcf);
        writeVideo(v,M)
    end
    close(v)
    
    v = VideoWriter([save_name simType ' KDE ' DATE ' rho_2 1D']);
    v.FrameRate = 5;
    open(v)
    
    figure
    for i=1:size(data1D,2)
        plot(X(1,:),data1D{2,i})
        axis([0 Lx 0 max(union(data1D{2,1},data1D{2,end}),[],'all')])
        xlabel('x')
        ylabel('Density')
        M = getframe(gcf);
        writeVideo(v,M)
    end
    close(v)

% %% Combined version
% 
%     PDF = cell(1,length(data));
%     PDF_1d = cell(1,length(data));
%     
%     p = Progress(length(data)) % Starts progress bar in Command Window.
%     for i = 1:length(data)
%         p.d(i) % Progress bar update
% %         i
%         Pos = data{i}.pos;
%         
%         [~,PDF{i},X,Y] = kde2d(Pos,N,[0 0],[Ly Lx],bwd(1),bwd(2));
%         PDF{i} = PDF{i}*N^2;
%         
%         PDF_1d{i} = trapz(Y(:,1),PDF{i},1)/Ly;
% 
% 
%     end
%     p.done() % Carrage return for progress bar display in Command Window.
%     
%     save('Shockwave Test after kde2d (07-05-2021) combined','PDF','PDF_1d','X','Y')
% 
% 
%     %% 2D movie
%     v = VideoWriter('Shockwave Test 2D after kde2d (07-05-2021) combined');
%     v.FrameRate = 5;
%     open(v)
%     
%     figure
%     for i=1:size(PDF,2)
%         mesh(X,Y,PDF{i})
%         axis([0 Lx 0 Ly 0 max(PDF{1},[],'all')])
%         xlabel('x')
%         ylabel('y')
%         zlabel('Density')
%         M = getframe(gcf);
%         writeVideo(v,M)
%     end
%     close(v)
%     
%     
%     %% 1D movie
%     v = VideoWriter('Shockwave Test 1D after kde2d (07-05-2021) combined');
%     v.FrameRate = 5;
%     open(v)
%     
%     figure
%     for i=1:size(PDF_1d,2)
%         plot(X(1,:),PDF_1d{i})
%         axis([0 Lx 0 max(PDF_1d{1},[],'all')])
%         xlabel('x')
%         ylabel('Density')
%         M = getframe(gcf);
%         writeVideo(v,M)
%     end
%     close(v)
    


end
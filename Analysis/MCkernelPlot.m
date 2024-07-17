function [PDF,PDF2,X,Y] = MCkernelPlot(data,params,varargin)
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
    addOptional(p,'Debug', 0, @islogical);
    addOptional(p,'StandardDevX', sqrt(0.007), @isnumeric); % for using exp func %sqrt(0.005)
%     addOptional(p,'StandardDevX', 0.005, @isnumeric); % for using mvnpdf
    addOptional(p,'StandardDevO', 1, @isnumeric);
    addOptional(p,'LengthX', 1, @isnumeric);
    addOptional(p,'LengthY',1,@isnumeric);
    addRequired(p,'MCdata');
    addRequired(p,'MCparams');
    
    parse(p,data,params, varargin{:} );
    
    
    Lx = params.grid*params.lengthG;
    Ly = Lx;
    dx = 0.01;
    dy = 0.01;
%     dt = params.dt;
    
%     sigmaX = p.Results.StandardDevX;
    sigmaX = sigmaCalc(params.N);
%     sigmaX = sigma;
    
    x = linspace(0,Lx,ceil(Lx/dx));
    y = linspace(0,Ly,ceil(Ly/dy));
    [X,Y] = meshgrid(x,y);
    xVec = [X(:) Y(:)];
    
%     Ratio = cell(1,params.Nt-2);
%     Diff = cell(1,params.Nt-2);
    PDF = cell(1,params.Nt);
%     xDeriv = cell(1,params.Nt);
%     yDeriv = cell(1,params.Nt);
%     tDeriv = cell(1,params.Nt);
    
    figure
    f1 = gcf;
    
    p = Progress(params.Nt) % Starts progress bar in Command Window.
%     for i = 1:20:params.Nt
    for i = 1:params.Nt
        p.d(i) % Progress bar update
%         i
        Pos = data{i}.pos;
        Ang = data{i}.o;
        tCurr = data{i}.t;
        
        [PDF{i}] = PDFcalc(xVec,sigmaX,Pos,Ang,Lx,Ly,params,tCurr,f1);
%         PDF{i} = reshape(LHS{i},length(y),length(x));
        
        PDFt = PDF{i};
%         PDF2{1,ceil(i/20)} = reshape(PDFt(1,:),sqrt(length(PDFt(1,:))),sqrt(length(PDFt(1,:))));
%         PDF2{2,ceil(i/20)} = reshape(PDFt(2,:),sqrt(length(PDFt(2,:))),sqrt(length(PDFt(2,:))));
        PDF2{1,i} = reshape(PDFt(1,:),sqrt(length(PDFt(1,:))),sqrt(length(PDFt(1,:))));
        PDF2{2,i} = reshape(PDFt(2,:),sqrt(length(PDFt(2,:))),sqrt(length(PDFt(2,:))));

    end
    p.done() % Carrage return for progress bar display in Command Window.
    
    save('Shockwave Test 1 after kernel density (06-10-2021)','PDF','PDF2')
    
%     v = VideoWriter('2020 MC sim (02-15-2021) rho_1 video 33pi-2, 66pi-2');
%     v.FrameRate = 5;
%     open(v)
%     
%     figure
%     for i=1:size(PDF2,2)
%         mesh(X,Y,PDF2{1,i})
%         axis([0 1 0 1 0 1])
%         M = getframe(gcf);
%         writeVideo(v,M)
%     end
%     close(v)
%     
%     v = VideoWriter('2020 MC sim (02-15-2021) rho_1 video 33pi-2, 66pi-2');
%     v.FrameRate = 5;
%     open(v)
%     
%     figure
%     for i=1:size(PDF2,2)
%         mesh(X,Y,PDF2{2,i})
%         axis([0 1 0 1 0 1])
%         M = getframe(gcf);
%         writeVideo(v,M)
%     end
%     close(v)


end

function [PDF] = PDFcalc(xVec,sigmaX,Pos,Ang,Lx,Ly,params,tCurr,f1)
    
    angles = sort(unique(Ang));
    PDF = zeros(length(angles),length(xVec));    
%     xDeriv = zeros(length(angles),length(xVec));
%     yDeriv = zeros(length(angles),length(xVec));
    
%     sigma = sigmaX*eye(2);

    for i = 1:length(angles)
%         i
        angCurr = Ang(Ang==angles(i));
        PosCurr = Pos(Ang==angles(i),:);
%         sigmaX = sigmaCalc(length(angCurr));
        for j = 1:length(angCurr)
            mu = PosCurr(j,:);
            for ii = -1:1
                for jj = -1:1
                    L = [Lx*ii Ly*jj];
%                     add = mvnpdf(xVec+repmat(L,length(xVec),1),mu,sigma)';
                    add = 1/(2*pi*sigmaX^2)*exp(-((xVec(:,1)'+L(1)-mu(1)).^2+(xVec(:,2)'+L(2)-mu(2)).^2)/(2*sigmaX^2));
                    PDF(i,:) = PDF(i,:) + add;
%                     xDeriv(i,:) = xDeriv(i,:)-(xVec(:,1)'+L(1)-mu(1))./sigmaX^2.*add;
%                     yDeriv(i,:) = yDeriv(i,:)-(xVec(:,2)'+L(2)-mu(2))./sigmaX^2.*add;
                end
            end
        end
        PDF(i,:) = PDF(i,:)/params.N;
%         xDeriv(i,:) = xDeriv(i,:)/params.N;
%         yDeriv(i,:) = yDeriv(i,:)/params.N;
        
        Z = reshape(PDF(i,:),sqrt(length(PDF(i,:))),sqrt(length(PDF(i,:))));
%         Zx = reshape(xDeriv(i,:),sqrt(length(xDeriv(i,:))),sqrt(length(xDeriv(i,:))));
%         Zy = reshape(yDeriv(i,:),sqrt(length(yDeriv(i,:))),sqrt(length(yDeriv(i,:))));
        
%         Z = boxAverage(Z);
%         Zx = boxAverage(Zx);
%         Zy = boxAverage(Zy);
%         
%         PDF(i,:) = reshape(Z,1,length(PDF(i,:)));
%         xDeriv(i,:) = reshape(Zx,1,length(xDeriv(i,:)));
%         yDeriv(i,:) = reshape(Zy,1,length(yDeriv(i,:)));
        
        cla(f1)
        X = reshape(xVec(:,1),sqrt(length(xVec(:,1))),sqrt(length(xVec(:,1))));
        Y = reshape(xVec(:,2),sqrt(length(xVec(:,2))),sqrt(length(xVec(:,2))));
        
        mesh(X,Y,Z)
        axis([0 1 0 1 0 5]);
        xlabel('x')
        ylabel('y')
        title(num2str(angles(i)))
        drawnow
        pause(1/2)
                
%         mesh(X,Y,Zx)
%         axis([0 1 0 1 -6 6]);
%         xlabel('x')
%         ylabel('y')
%         title(strcat(num2str(angles(i)),' x derivative'))
%         drawnow
%         pause(1/2)
%         
%         mesh(X,Y,Zy)
%         axis([0 1 0 1 -5 5]);
%         xlabel('x')
%         ylabel('y')
%         title(strcat(num2str(angles(i)),' y derivative'))
%         drawnow
%         pause(1/2)
        
    end
    
    
end


function sigma = sigmaCalc(N)
    
    initSigma = 0.03; % Taken from in MC code
    
%     sigma = sqrt(initSigma)*(32/(3*sqrt(2)*N))^(1/6);

    sigma = sqrt(initSigma)*N^(-1/6);
    
end
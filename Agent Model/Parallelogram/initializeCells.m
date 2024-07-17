 function [Pos,Ang,BinLocations] = initializeCells(params,IC_type,varargin)

% Initializes cell positions and angle based on the specified initial
% condition type.

% Inputs:
%
%   params - A structure containing the parameters used in  the main
%   simulation.
%
%   IC_type - A string indicating which initialization type to run.
%
% Outputs:
%
%   Pos - The spatial locations of the initialized cell as a N x 2 array.
%
%   Ang - The orientations of the initialized cells, each given as an angle
%   between -pi and pi.
%
%   BinLocations - The locations in the interval (-pi,pi) where bins edges
%   are defined. Used in color plotting in main script.

    p = inputParser;
    addOptional(p,'InitialMean',[],@isnumeric);
    addOptional(p,'InitialSD',[],@isnumeric)
    addOptional(p,'InitialProp',[],@isnumeric);
    parse(p,varargin{:});  

    
    if strcmp(IC_type,'Uniform')
        % Uniform initial conditions
        prop1 = 0.8; %Proportion of cells in main cone
        prop2 = 1-prop1; %Proportion of cells in secondary cone
        rho1 = 2/3; %Proportion of cells in main cone moving with Ang1
        rho2 = 1/3; %Analogue of rho1 for secondary cone
        xinit = params.grid*params.lengthGx*rand(params.N,1);
        yinit = params.grid*params.lengthGy*rand(params.N,1);
        Pos = [xinit,yinit];
        
        % Impose periodic boundary conditions (just in case)
        Pos = mod(Pos,params.grid*[params.lengthGx,params.lengthGy]);
        
        % Uniform over small angle interval
        Ang1 = pi/3*rand(floor(prop1*rho1*params.N),1) - pi/6;
        Ang2 = pi/3*rand(ceil(params.N*prop1)-floor(prop1*rho1*params.N),1)+5*pi/6;
        Ang3 = 2*pi/3*rand(floor(prop2*rho2*params.N),1)+pi/6;
        Ang4 = 2*pi/3*rand(ceil(prop2*params.N)-floor(prop2*rho2*params.N),1)+7*pi/6;
        Ang = cat(1,Ang1,Ang2,Ang3,Ang4)-5*pi/6;
        
        BinLocations = [-pi/6,pi/6,5*pi/6,7*pi/6,11*pi/6]-5*pi/6;

        % Gaussian setup
    %     rho = 1/2; %Controls proportion of cells in the two gaussians
    %     mu1 = [0.25 0.5];
    %     sigma1 = 0.01*eye(2);
    %     mu2 = [0.75 0.5];
    %     sigma2 = 0.01*eye(2);

        % multivariate
    %     Pos1 = mvnrnd(mu1,sigma1,floor(rho*params.N));
    %     Pos2 = mvnrnd(mu2,sigma2,params.N-floor(rho*params.N));
    %     Pos = [Pos1;Pos2];

            % univariate in x-direction, uniform in y-direction
    %     Pos1 = mu1(1) + sigma1(1,1)*randn(floor(rho*params.N),1);
    %     Pos2 = mu2(1) + sigma2(1,1)*randn(params.N-floor(rho*params.N),1);
    %     Pos = [Pos1 rand(floor(rho*params.N),1)*params.grid*params.lengthG ...
    %         ;Pos2 rand(params.N-floor(rho*params.N),1)*params.grid*params.lengthG];

        

        % Plots initial conditions to check periodicity at boundaries.
            % multivariate
    %     Z1 = mvnpdf(Pos(1:length(Pos1),:),mu1,sigma1);
    %     Z2 = mvnpdf(Pos(length(Pos1)+1:end,:),mu2,sigma2);
            % univariate in x-direction, uniform in y-direction
    %     Z1 = normpdf(Pos1,mu1(1),sigma1(1,1));
    %     Z2 = normpdf(Pos2,mu2(1),sigma2(1,1));

    %     scatter3(Pos(1:length(Pos1),1),Pos(1:length(Pos1),2),Z1,'r')
    %     hold on
    %     scatter3(Pos(length(Pos1)+1:end,1),Pos(length(Pos1)+1:end,2),Z2,'b')
    %     pause(2)
    %     waitforbuttonpress

    %     xinit1 = params.grid*params.lengthG/2*rand(params.N/2,1);
    %     xinit2 = params.grid*params.lengthG/2*(1+rand(params.N/2,1));
    %     yinit = params.grid*params.lengthG*rand(params.N,1);
    %     Pos = [[xinit1;xinit2],yinit];

        % Initialize orientation vectors. Set to be uniformly distributed.
    %     Ang = 2*pi*rand(params.N,1)-pi;

        % Fixed vectors.
    %     Ang1 = -0.1*ones(floor(2*params.N/3),1);
    %     Ang2 = pi/4*ones(params.N-floor(2*params.N/3),1);
    %     Ang = [Ang1;Ang2];
    %     BinLocations = [-pi,0,pi];

        % Fixed vectors #2.
    %     Ang1 = 3*pi/4*ones(floor(params.N/2),1);
    %     Ang2 = 1*pi/4*ones(params.N-floor(params.N/2),1);
    %     Ang = [Ang1;Ang2];
    %     BinLocations = [0,pi/2,pi];

        % Fixed vectors #3
    %     Ang1 = 1*pi/6*ones(floor(rho*params.N),1);
    %     Ang2 = 5*pi/6*ones(params.N-floor(rho*params.N),1);
    %     Ang = [Ang1;Ang2];
%         BinLocations = [-pi/6,pi/6,5*pi/6,7*pi/6,11*pi/6]-5*pi/6;

    %     BinLocations = linspace(0,2*pi,6);


    elseif strcmp(IC_type,'UniformBand')
        if ~isempty(p.Results.InitialProp)
            rho = p.Results.InitialProp;
        else
            rho = 1/2; %Controls proportion of cells in the two uniform dists
        end
        xinit = params.grid*params.lengthGx*rand(params.N,1);
        yinit = params.grid*params.lengthGy*rand(params.N,1);
        Pos = [xinit,yinit];
        
        % Impose periodic boundary conditions
        Pos = mod(Pos,params.grid*[params.lengthGx,params.lengthGy]);

        Ang1 = 1*pi/4*ones(floor(params.N/2),1);
        Ang2 = 3*pi/4*ones(params.N-floor(params.N/2),1);
        Ang = [Ang1;Ang2];

    BinLocations = [0,pi/2,pi];
    
    elseif strcmp(IC_type,'UniformTrue')
        % Uniform initial conditions
        xinit = params.grid*params.lengthGx*rand(params.N,1);
        yinit = params.grid*params.lengthGy*rand(params.N,1);
        Pos = [xinit,yinit];
        
        % Impose periodic boundary conditions (just in case)
        Pos = mod(Pos,params.grid*[params.lengthGx,params.lengthGy]);
        
        % Uniform over small angle interval
        Ang = 2*pi*rand(params.N,1)-pi;
        
        BinLocations = [-pi,-pi/2,0,pi/2,pi];
        
    elseif strcmp(IC_type,'UniformBias')
        % Uniform initial conditions
        prop1 = 0.5; %Proportion of cells in main cone
        prop2 = 1-prop1; %Proportion of cells in secondary cone
        xinit = params.grid*params.lengthGx*rand(params.N,1);
        yinit = params.grid*params.lengthGy*rand(params.N,1);
        Pos = [xinit,yinit];
        
        % Impose periodic boundary conditions (just in case)
        Pos = mod(Pos,params.grid*[params.lengthGx,params.lengthGy]);
        
        % Gaussian over angle intervals
        Ang1 = pi/3*rand(floor(prop1*params.N),1);
        Ang2 = pi/3*rand(params.N-floor(prop1*params.N),1)-pi/2-pi/6;
        Ang = cat(1,Ang1,Ang2);
        
        BinLocations = [-pi,-pi/2,0,pi/2,pi];
        
    elseif strcmp(IC_type,'GaussianBias')
        % Uniform initial conditions
        prop1 = 0.3; %Proportion of cells in main cone
        prop2 = 1-prop1; %Proportion of cells in secondary cone
        xinit = params.grid*params.lengthGx*rand(params.N,1);
        yinit = params.grid*params.lengthGy*rand(params.N,1);
        Pos = [xinit,yinit];
        
        % Impose periodic boundary conditions (just in case)
        Pos = mod(Pos,params.grid*[params.lengthGx,params.lengthGy]);
        
        % Gaussian over angle intervals
        Ang1 = mod(pi/18*randn(floor(prop1*params.N),1)+pi/6+pi,2*pi)-pi;
        Ang2 = mod(pi/18*randn(params.N-floor(prop1*params.N),1)-pi/2+pi,2*pi)-pi;
        Ang = cat(1,Ang1,Ang2);
        
        BinLocations = [-pi,-pi/2,0,pi/2,pi];
    
    elseif strcmp(IC_type,'Test')
        % N = 2, K = 0.1, dt = 0.0025
%         Pos = [0.38 0.501;...
%                 0.62 0.5];
%         Ang = [0-0.01;pi];
%         Pos = [0.38 0.51;...
%                 0.62 0.5];
%         Ang = [0-0.1;pi];
%         Pos = [0.38 0.51;...
%                 0.62 0.5];
%         Ang = [0-0.05;pi];
%         BinLocations = [-pi/2, pi/2, 3*pi/2];

        Pos = [0.5-1/sqrt(3)*0.5*0.2 0.5;...
                0.5+1/sqrt(3)*0.5*0.2 0.5-0.001]*params.lengthGx;
        Ang = [pi/3;2*pi/3];
        BinLocations = [-pi/2, pi/2, 3*pi/2];

        Pos = [0.5-1/sqrt(2)*0.015*0.2 0.5;...
                0.5+1/sqrt(2)*0.015*0.2 0.5-0]*params.lengthGx;
        Ang = [pi/4;3*pi/4];
        BinLocations = [-pi/2, pi/2, 3*pi/2];

%         Pos = [0.5-1/sqrt(3)*0.5*0.2 0.5;...
%                 0.5+1/sqrt(3)*0.5*0.2 0.5;...
%                 0.5+1/sqrt(3)*0.5*0.2-0.05 0.5-0.05*tan(pi/3)];
%         Ang = [pi/3;2*pi/3;2*pi/3];
%         BinLocations = [-pi/2, pi/2, 3*pi/2];

    elseif strcmp(IC_type,'Lattice_1_R0')
        Nm = params.lengthGx;
        y = randi([1,params.lengthGy],params.N,1);
        x = zeros(params.N,1);
        
        % left-moving band
        for i=params.N/2+1:params.N
            x(i) = normrnd(0,floor(25/400*Nm));
            while x(i) < -floor(149/400*Nm) || x(i) > floor(350/400*Nm)
                x(i) = normrnd(0,floor(25/400*Nm));
            end
            x(i) = floor(x(i))+floor(150/400*Nm);
        end

        % right-moving band
        for i=1:params.N/2
            x(i) = normrnd(0,floor(25/400*Nm));
            while x(i) < -floor(249/400*Nm) || x(i) > floor(150/400*Nm)
                x(i) = normrnd(0,floor(25/400*Nm));
            end
            x(i) = floor(x(i))+floor(250/400*Nm);
        end

    Pos = [x,y];

    % Impose periodic boundary conditions
    Pos = mod(Pos,params.grid*[params.lengthGx,params.lengthGy]);

%     Ang1 = 3*pi/4*ones(floor(params.N/2),1);
%     Ang2 = 1*pi/4*ones(params.N-floor(params.N/2),1);
    Ang1 = (3*pi/4-pi/8)*ones(floor(params.N/2),1);
    Ang2 = (pi/4+pi/8)*ones(params.N-floor(params.N/2),1);
    Ang = [Ang1;Ang2];

    BinLocations = [0,pi/2,pi];
    
    
    elseif strcmp(IC_type,'Lattice_1_R1')
        Nr = 1;
        Nm = params.lengthGx;
        y = randi([1,params.lengthGy*2^Nr],params.N,1)/2^Nr;
        x = zeros(params.N,1);
        
        % left-moving band
        for i=params.N/2+1:params.N
            x(i) = normrnd(0,floor(25/400*Nm*2^Nr));
            while x(i) < -floor(149/400*Nm*2^Nr) || x(i) > floor(350/400*Nm*2^Nr)
                x(i) = normrnd(0,floor(25/400*Nm*2^Nr));
            end
            x(i) = floor(x(i))/2^Nr+floor(150/400*Nm);
        end

        % right-moving band
        for i=1:params.N/2
            x(i) = normrnd(0,floor(25/400*Nm*2^Nr));
            while x(i) < -floor(249/400*Nm*2^Nr) || x(i) > floor(150/400*Nm*2^Nr)
                x(i) = normrnd(0,floor(25/400*Nm*2^Nr));
            end
            x(i) = floor(x(i))/2^Nr+floor(250/400*Nm);
        end

    Pos = [x,y];

    % Impose periodic boundary conditions
    Pos = mod(Pos,params.grid*[params.lengthGx,params.lengthGy]);

    Ang1 = 3*pi/4*ones(floor(params.N/2),1);
    Ang2 = 1*pi/4*ones(params.N-floor(params.N/2),1);
    Ang = [Ang1;Ang2];

    BinLocations = [0,pi/2,pi];
    
    
    elseif strcmp(IC_type,'Lattice_1_R2')
        Nr = 2;
        Nm = params.lengthGx;
        y = randi([1,params.lengthGy*2^Nr],params.N,1)/2^Nr;
        x = zeros(params.N,1);
        
        % left-moving band
        for i=params.N/2+1:params.N
            x(i) = normrnd(0,floor(25/400*Nm*2^Nr));
            while x(i) < -floor(149/400*Nm*2^Nr) || x(i) > floor(350/400*Nm*2^Nr)
                x(i) = normrnd(0,floor(25/400*Nm*2^Nr));
            end
            x(i) = floor(x(i))/2^Nr+floor(150/400*Nm);
        end

        % right-moving band
        for i=1:params.N/2
            x(i) = normrnd(0,floor(25/400*Nm*2^Nr));
            while x(i) < -floor(249/400*Nm*2^Nr) || x(i) > floor(150/400*Nm*2^Nr)
                x(i) = normrnd(0,floor(25/400*Nm*2^Nr));
            end
            x(i) = floor(x(i))/2^Nr+floor(250/400*Nm);
        end

    Pos = [x,y];

    % Impose periodic boundary conditions
    Pos = mod(Pos,params.grid*[params.lengthGx,params.lengthGy]);

    Ang1 = 3*pi/4*ones(floor(params.N/2),1);
    Ang2 = 1*pi/4*ones(params.N-floor(params.N/2),1);
    Ang = [Ang1;Ang2];

    BinLocations = [0,pi/2,pi];
    
    
    
    
    elseif strcmp(IC_type,'Lattice_2_R0')
        Nm = params.lengthGx;
        y = randi([1,params.lengthGy],params.N,1);
        x = zeros(params.N,1);
        
        Pr1 =zeros(1,Nm); % CDF for first "square" band
        Pr2 =zeros(1,Nm); % CDF for second "square" band
        f1 = @(x) exp(-((x - floor(150/400*Nm))/floor(50/400*Nm)).^8);
        f2 = @(x) exp(-((x - floor(250/400*Nm))/floor(50/400*Nm)).^8);
        for i=1:Nm
            Pr1(i) = integral(f1,0,i);
            Pr2(i) = integral(f2,0,i);
        end
        Pr1 = Pr1/integral(f1,0,Nm);
        Pr2 = Pr2/integral(f2,0,Nm);
        
        % left-moving band
        for i=params.N/2+1:params.N
            xtemp = rand;
            k = 1;
            while xtemp > Pr1(k)
                k = k+1;
            end
            x(i) = k;
        end

        % right-moving band
        for i=1:params.N/2
            xtemp = rand;
            k = 1;
            while xtemp > Pr2(k)
                k = k+1;
            end
            x(i) = k;
        end

    Pos = [x,y];

    % Impose periodic boundary conditions
    Pos = mod(Pos,params.grid*[params.lengthGx,params.lengthGy]);

    Ang1 = 3*pi/4*ones(floor(params.N/2),1);
    Ang2 = 1*pi/4*ones(params.N-floor(params.N/2),1);
    Ang = [Ang1;Ang2];

    BinLocations = [0,pi/2,pi];
    
    
    elseif strcmp(IC_type,'Lattice_2_R1')
        Nr = 1;
        Nm = params.lengthGx;
        y = randi([1,params.lengthGy]*2^Nr,params.N,1)/2^Nr;
        x = zeros(params.N,1);
        
        Pr1 =zeros(1,Nm*2^Nr); % CDF for first "square" band
        Pr2 =zeros(1,Nm*2^Nr); % CDF for second "square" band
        f1 = @(x) exp(-((x - floor(150/400*Nm*2^Nr))/floor(50/400*Nm*2^Nr)).^8);
        f2 = @(x) exp(-((x - floor(250/400*Nm*2^Nr))/floor(50/400*Nm*2^Nr)).^8);
        for i=1:Nm*2^Nr
            Pr1(i) = integral(f1,0,i);
            Pr2(i) = integral(f2,0,i);
        end
        Pr1 = Pr1/integral(f1,0,Nm*2^Nr);
        Pr2 = Pr2/integral(f2,0,Nm*2^Nr);
        
        % left-moving band
        for i=params.N/2+1:params.N
            xtemp = rand;
            k = 1;
            while xtemp > Pr1(k)
                k = k+1;
            end
            x(i) = k/2^Nr;
        end

        % right-moving band
        for i=1:params.N/2
            xtemp = rand;
            k = 1;
            while xtemp > Pr2(k)
                k = k+1;
            end
            x(i) = k/2^Nr;
        end

    Pos = [x,y];

    % Impose periodic boundary conditions
    Pos = mod(Pos,params.grid*[params.lengthGx,params.lengthGy]);

    Ang1 = 3*pi/4*ones(floor(params.N/2),1);
    Ang2 = 1*pi/4*ones(params.N-floor(params.N/2),1);
    Ang = [Ang1;Ang2];

    BinLocations = [0,pi/2,pi];
    
    
    elseif strcmp(IC_type,'Lattice_2_R2')
        Nr = 2;
        Nm = params.lengthGx;
        y = randi([1,params.lengthGy]*2^Nr,params.N,1)/2^Nr;
        x = zeros(params.N,1);
        
        Pr1 =zeros(1,Nm*2^Nr); % CDF for first "square" band
        Pr2 =zeros(1,Nm*2^Nr); % CDF for second "square" band
        f1 = @(x) exp(-((x - floor(150/400*Nm*2^Nr))/floor(50/400*Nm*2^Nr)).^8);
        f2 = @(x) exp(-((x - floor(250/400*Nm*2^Nr))/floor(50/400*Nm*2^Nr)).^8);
        for i=1:Nm*2^Nr
            Pr1(i) = integral(f1,0,i);
            Pr2(i) = integral(f2,0,i);
        end
        Pr1 = Pr1/integral(f1,0,Nm*2^Nr);
        Pr2 = Pr2/integral(f2,0,Nm*2^Nr);
        
        % left-moving band
        for i=params.N/2+1:params.N
            xtemp = rand;
            k = 1;
            while xtemp > Pr1(k)
                k = k+1;
            end
            x(i) = k/2^Nr;
        end

        % right-moving band
        for i=1:params.N/2
            xtemp = rand;
            k = 1;
            while xtemp > Pr2(k)
                k = k+1;
            end
            x(i) = k/2^Nr;
        end

    Pos = [x,y];

    % Impose periodic boundary conditions
    Pos = mod(Pos,params.grid*[params.lengthGx,params.lengthGy]);

    Ang1 = 3*pi/4*ones(floor(params.N/2),1);
    Ang2 = 1*pi/4*ones(params.N-floor(params.N/2),1);
    Ang = [Ang1;Ang2];

    BinLocations = [0,pi/2,pi];
    
    
    
    elseif strcmp(IC_type,'Flat_Gaussian')
        Nr = 3;
        Nm = params.lengthGx;
        y = params.lengthGy*rand(params.N,1);
        x = zeros(params.N,1);
        
        Pr1 =zeros(1,Nm*2^Nr); % CDF for first "square" band
        Pr2 =zeros(1,Nm*2^Nr); % CDF for second "square" band
        f1 = @(x) exp(-((x - floor(p.Results.InitialMean(1)*Nm*2^Nr))/floor(2*p.Results.InitialSD(1)*Nm*2^Nr)).^8);
        f2 = @(x) exp(-((x - floor(p.Results.InitialMean(2)*Nm*2^Nr))/floor(2*p.Results.InitialSD(2)*Nm*2^Nr)).^8);
        for i=1:Nm*2^Nr
            Pr1(i) = integral(f1,0,i);
            Pr2(i) = integral(f2,0,i);
        end
        Pr1 = Pr1/integral(f1,0,Nm*2^Nr);
        Pr2 = Pr2/integral(f2,0,Nm*2^Nr);
        
        % left-moving band
        for i=params.N/2+1:params.N
            xtemp = rand;
            k = 1;
            while xtemp > Pr1(k)
                k = k+1;
            end
            x(i) = (k+(rand-0.5))/2^Nr;
        end

        % right-moving band
        for i=1:params.N/2
            xtemp = rand;
            k = 1;
            while xtemp > Pr2(k)
                k = k+1;
            end
            x(i) = (k+(rand-0.5))/2^Nr;
        end

    Pos = [x,y];

    % Impose periodic boundary conditions
    Pos = mod(Pos,params.grid*[params.lengthGx,params.lengthGy]);

    Ang1 = 3*pi/4*ones(floor(params.N/2),1);
    Ang2 = 1*pi/4*ones(params.N-floor(params.N/2),1);
    Ang = [Ang1;Ang2];

    BinLocations = [0,pi/2,pi];
    
    
        
        
    elseif strcmp(IC_type,'Shockwave')
        % alignMCwavetest.m
        % May Wave Test initial conditions
        prop1 = 1/4; %Proportion of cells in for rho1 across shock L to R
        prop2 = 3/5; %Proportion of cells in for rho2 across shock L to R
        propRho12 = 2/(2+5/3); %Proportion of cells int(rho1)/[int(rho1)+int(rho2)] in Left;
        rho1 = 2/3; %Proportion of cells in main cone moving with Ang1
        rho2 = 1/3; %Analogue of rho1 for secondary code
        xinit1L = params.grid*params.lengthGx/2*rand(floor(prop1*propRho12*params.N),1);
        xinit1R = params.grid*params.lengthGx/2*rand(ceil(propRho12*params.N)...
            -floor(prop1*propRho12*params.N),1)...
            + params.grid*params.lengthGx/2;
        xinit2L = params.grid*params.lengthGx/2*rand(floor(prop2*(1-propRho12)*params.N),1);
        xinit2R = params.grid*params.lengthGx/2*rand(ceil((1-propRho12)*params.N)...
            -floor(prop2*(1-propRho12)*params.N),1)...
            + params.grid*params.lengthGx/2;
        extraCells = length(xinit1L)+length(xinit1R)+length(xinit2L)+length(xinit2R) - params.N;
        if extraCells > 0
            for i = 1:extraCells
                if mod(i,4) == 0
                    xinit1L = xinit1L(1:end-1);
                elseif mod(i,4) == 1
                    xinit1R = xinit1R(1:end-1);
                elseif mod(i,4) == 2
                    xinit2L = xinit2L(1:end-1);
                elseif mod(i,4) == 3
                    xinit2R = xinit2R(1:end-1);
                end
            end
        end
        yinit = params.grid*params.lengthGy*rand(params.N,1);

        Pos = [[xinit1L;xinit1R;xinit2L;xinit2R],yinit];
        
        % Fixed vectors #2.
        Ang1 = 3*pi/4*ones(floor(params.N/2),1);
        Ang2 = 1*pi/4*ones(params.N-floor(params.N/2),1);
        Ang = [Ang1;Ang2];
        
        BinLocations = [0,pi/2,pi];


    elseif strcmp(IC_type,'Shockwave2')
        % alignMCwavetest.m
        % May Wave Test initial conditions
        prop1 = 1/(1+6); %Proportion of cells in for rho1 across shock L to R
        prop2 = 4/(4+1.887); %Proportion of cells in for rho2 across shock L to R
        propRho12 = 3.5/(3.5+2+1.887/2); %Proportion of cells int(rho1)/[int(rho1)+int(rho2)] in Left;
        xinit1L = params.grid*params.lengthGx/2*rand(floor(prop1*propRho12*params.N),1);
        xinit1R = params.grid*params.lengthGx/2*rand(ceil(propRho12*params.N)...
            -floor(prop1*propRho12*params.N),1)...
            + params.grid*params.lengthGx/2;
        xinit2L = params.grid*params.lengthGx/2*rand(floor(prop2*(1-propRho12)*params.N),1);
        xinit2R = params.grid*params.lengthGx/2*rand(ceil((1-propRho12)*params.N)...
            -floor(prop2*(1-propRho12)*params.N),1)...
            + params.grid*params.lengthGx/2;
        extraCells = length(xinit1L)+length(xinit1R)+length(xinit2L)+length(xinit2R) - params.N;
        if extraCells > 0
            for i = 1:extraCells
                if mod(i,4) == 0
                    xinit1L = xinit1L(1:end-1);
                elseif mod(i,4) == 1
                    xinit1R = xinit1R(1:end-1);
                elseif mod(i,4) == 2
                    xinit2L = xinit2L(1:end-1);
                elseif mod(i,4) == 3
                    xinit2R = xinit2R(1:end-1);
                end
            end
        end
        yinit = params.grid*params.lengthGy*rand(params.N,1);

        Pos = [[xinit1L;xinit1R;xinit2L;xinit2R],yinit];
        
        % Fixed vectors #2.
        Ang1 = 3*pi/4*ones(floor(params.N/2),1);
        Ang2 = 1*pi/4*ones(params.N-floor(params.N/2),1);
        Ang = [Ang1;Ang2];
        
        BinLocations = [0,pi/2,pi];
        
    elseif strcmp(IC_type,'Squarewave')
        % Uniform wave setup
        if ~isempty(p.Results.InitialProp)
            rho = p.Results.InitialProp;
        else
            rho = 1/2; %Controls proportion of cells in the two gaussians
        end
        if ~isempty(p.Results.InitialSD)
            sigma1 = p.Results.InitialSD(1); % Base width of rectangle 1
            mu1 = 0.5-sigma1; % Location of left base of rectangle 1
            sigma2 = p.Results.InitialSD(2); % Base width of rectangle 2
            mu2 = 0.5; % Location of left base of triangle 2
        else
            sigma1 = 0.25; % Base width of rectangle 1
            mu1 = 0.5-sigma1; % Location of left base of rectangle 1
            sigma2 = sigma1; % Base width of rectangle 2
            mu2 = 0.5; % Location of left base of rectangle 2
        end

            % univariate in x-direction, uniform in y-direction
        xinit1 = params.grid*params.lengthGx*(mu1 + sigma1*(rand(floor(rho*params.N),1)));
        xinit2 = params.grid*params.lengthGx*(mu2 + sigma2*(rand(params.N-floor(rho*params.N),1)));
        yinit = params.grid*params.lengthGy*rand(params.N,1);
        Pos = [[xinit1;xinit2],yinit];
        
        % Impose periodic boundary conditions
        Pos = mod(Pos,params.grid*[params.lengthGx,params.lengthGy]);
        
        % Fixed vectors #2.
        Ang1 = 1*pi/4*ones(floor(params.N/2),1);
        Ang2 = 3*pi/4*ones(params.N-floor(params.N/2),1);
        Ang = [Ang1;Ang2];
        
        BinLocations = [0,pi/2,pi];
        

    elseif strcmp(IC_type,'AsymWave')
            % Asym Gaussian wave setup
%         mu1 = [0.30 0.5];
%         sigma1 = 0.05*eye(2);
%         mu2 = [0.70 0.5];
%         sigma2 = 0.05*eye(2);
        if ~isempty(p.Results.InitialProp)
            rho = p.Results.InitialProp;
            back_rho = p.Results.BackProp;
        else
            gauss_rho = 1/15; %Controls proportion of cells in the gaussians
            back_rho = 7/15;
        end
        if ~isempty(p.Results.InitialSD) && isempty(p.Results.InitialMean)
            sigma1 = p.Results.InitialSD(1)*params.lengthGx; % SD of left gaussian
            mu1 = 0.5*params.lengthGx-5*sigma1; % Location of left gaussian
%             sigma2 = p.Results.InitialSD(2)*params.lengthGx; % SD of right gaussian
%             mu2 = 0.5*params.lengthGx+5*sigma2*params.lengthGx; % Location of right gaussian
        elseif ~isempty(p.Results.InitialSD) && ~isempty(p.Results.InitialMean)
            sigma1 = p.Results.InitialSD(1)*params.lengthGx; % SD of left gaussian
            mu1 = p.Results.InitialMean(1)*params.lengthGx; % Location of left gaussian
%             sigma2 = p.Results.InitialSD(2)*params.lengthGx; % SD of right gaussian
%             mu2 = p.Results.InitialMean(2)*params.lengthGx; % Location of right gaussian
        else
            mu1 = [0.32 0.5]*params.lengthGx;
            sigma1 = 0.05*eye(2)*params.lengthGx;
%             mu2 = [0.68 0.5]*params.lengthGx;
%             sigma2 = 0.05*eye(2)*params.lengthGx;
        end


            % univariate in x-direction, uniform in y-direction
        Pos1 = mu1(1) + sigma1(1,1)*randn(round(gauss_rho*params.N),1);
%         Pos2 = mu2(1) + sigma2(1,1)*randn(round((1-rho)*params.N*(1-back_rho)),1);
        Pos1b = params.lengthGx*rand(round(back_rho*params.N),1);
        Pos2b = params.lengthGx*rand(round(back_rho*params.N),1);
        Posx = [Pos1; Pos1b; Pos2b];
        if length(Posx) >= params.N
            Pos = [Posx(1:params.N) rand(params.N,1)*params.lengthGy];
        elseif length(Posx) < params.N
            a = params.N - length(Posx);
            Pos = [[Posx; mu1(1) + sigma1(1,1)*randn(a,1)]...
                rand(params.N,1)*params.lengthGy];
        end
        
        % Impose periodic boundary conditions
        Pos = mod(Pos,params.grid*[params.lengthGx,params.lengthGy]);
        
        % Fixed vectors #2.
        Ang1 = 1*pi/4*ones(floor(params.N*(back_rho+gauss_rho)),1);
        Ang2 = 3*pi/4*ones(ceil(params.N*back_rho),1);
        Ang = [Ang1;Ang2];
        
        BinLocations = [0,pi/2,pi];

    elseif strcmp(IC_type,'RaisedWave')
            % Gaussian wave setup
%         mu1 = [0.30 0.5];
%         sigma1 = 0.05*eye(2);
%         mu2 = [0.70 0.5];
%         sigma2 = 0.05*eye(2);
        if ~isempty(p.Results.InitialProp)
            rho = p.Results.InitialProp;
            back_rho = p.Results.BackProp;
        else
            rho = 1/2; %Controls proportion of cells in the two gaussians
            back_rho = 7/8;
        end
        if ~isempty(p.Results.InitialSD) && isempty(p.Results.InitialMean)
            sigma1 = p.Results.InitialSD(1)*params.lengthGx; % SD of left gaussian
            mu1 = 0.5*params.lengthGx-5*sigma1; % Location of left gaussian
            sigma2 = p.Results.InitialSD(2)*params.lengthGx; % SD of right gaussian
            mu2 = 0.5*params.lengthGx+5*sigma2*params.lengthGx; % Location of right gaussian
        elseif ~isempty(p.Results.InitialSD) && ~isempty(p.Results.InitialMean)
            sigma1 = p.Results.InitialSD(1)*params.lengthGx; % SD of left gaussian
            mu1 = p.Results.InitialMean(1)*params.lengthGx; % Location of left gaussian
            sigma2 = p.Results.InitialSD(2)*params.lengthGx; % SD of right gaussian
            mu2 = p.Results.InitialMean(2)*params.lengthGx; % Location of right gaussian
        else
            mu1 = [0.32 0.5]*params.lengthGx;
            sigma1 = 0.05*eye(2)*params.lengthGx;
            mu2 = [0.68 0.5]*params.lengthGx;
            sigma2 = 0.05*eye(2)*params.lengthGx;
        end


            % univariate in x-direction, uniform in y-direction
        Pos1 = mu1(1) + sigma1(1,1)*randn(round(rho*params.N*(1-back_rho)),1);
        Pos2 = mu2(1) + sigma2(1,1)*randn(round((1-rho)*params.N*(1-back_rho)),1);
        Pos1b = params.lengthGx*rand(round(rho*params.N*back_rho),1);
        Pos2b = params.lengthGx*rand(round((1-rho)*params.N*back_rho),1);
        Posx = [Pos1; Pos1b; Pos2; Pos2b];
        if length(Posx) >= params.N
            Pos = [Posx(1:params.N) rand(params.N,1)*params.lengthGy];
        elseif length(Posx) < params.N
            a = params.N - length(Posx);
            Pos = [[Posx; mu1(1) + sigma1(1,1)*randn(a,1)]...
                rand(params.N,1)*params.lengthGy];
        end
        
        % Impose periodic boundary conditions
        Pos = mod(Pos,params.grid*[params.lengthGx,params.lengthGy]);
        
        % Fixed vectors #2.
        Ang1 = 1*pi/4*ones(floor(params.N/2),1);
        Ang2 = 3*pi/4*ones(params.N-floor(params.N/2),1);
        Ang = [Ang1;Ang2];
        
        BinLocations = [0,pi/2,pi];


    elseif strcmp(IC_type,'Wave')
            % Gaussian wave setup
%         mu1 = [0.30 0.5];
%         sigma1 = 0.05*eye(2);
%         mu2 = [0.70 0.5];
%         sigma2 = 0.05*eye(2);
        if ~isempty(p.Results.InitialProp)
            rho = p.Results.InitialProp;
        else
            rho = 1/2; %Controls proportion of cells in the two gaussians
        end
        if ~isempty(p.Results.InitialSD) && isempty(p.Results.InitialMean)
            sigma1 = p.Results.InitialSD(1)*params.lengthGx; % SD of left gaussian
            mu1 = 0.5*params.lengthGx-5*sigma1; % Location of left gaussian
            sigma2 = p.Results.InitialSD(2)*params.lengthGx; % SD of right gaussian
            mu2 = 0.5*params.lengthGx+5*sigma2*params.lengthGx; % Location of right gaussian
        elseif ~isempty(p.Results.InitialSD) && ~isempty(p.Results.InitialMean)
            sigma1 = p.Results.InitialSD(1)*params.lengthGx; % SD of left gaussian
            mu1 = p.Results.InitialMean(1)*params.lengthGx; % Location of left gaussian
            sigma2 = p.Results.InitialSD(2)*params.lengthGx; % SD of right gaussian
            mu2 = p.Results.InitialMean(2)*params.lengthGx; % Location of right gaussian
        else
            mu1 = [0.32 0.5]*params.lengthGx;
            sigma1 = 0.05*eye(2)*params.lengthGx;
            mu2 = [0.68 0.5]*params.lengthGx;
            sigma2 = 0.05*eye(2)*params.lengthGx;
        end


            % univariate in x-direction, uniform in y-direction
        Pos1 = mu1(1) + sigma1(1,1)*randn(floor(rho*params.N),1);
        Pos2 = mu2(1) + sigma2(1,1)*randn(params.N-floor(rho*params.N),1);
        Pos = [Pos1 rand(floor(rho*params.N),1)*params.grid*params.lengthGy ...
            ;Pos2 rand(params.N-floor(rho*params.N),1)*params.grid*params.lengthGy];
        
        % Impose periodic boundary conditions
        Pos = mod(Pos,params.grid*[params.lengthGx,params.lengthGy]);
        
        % Fixed vectors #2.
        Ang1 = 1*pi/4*ones(floor(params.N/2),1);
        Ang2 = 3*pi/4*ones(params.N-floor(params.N/2),1);
        Ang = [Ang1;Ang2];
        
        BinLocations = [0,pi/2,pi];
        
    elseif strcmp(IC_type,'Wedge')
        % Triangle wave setup
        if ~isempty(p.Results.InitialProp)
            rho = p.Results.InitialProp;
        else
            rho = 1/2; %Controls proportion of cells in the two gaussians
        end
        if ~isempty(p.Results.InitialSD)
            sigma1 = p.Results.InitialSD(1); % Base width of triangle 1
            mu1 = 0.5-sigma1; % Location of left base of triangle 1
            sigma2 = p.Results.InitialSD(2); % Base width of triangle 2
            mu2 = 0.5; % Location of left base of triangle 2
        else
            sigma1 = 0.25; % Base width of triangle 1
            mu1 = 0.5-sigma1; % Location of left base of triangle 1
            sigma2 = sigma1; % Base width of triangle 2
            mu2 = 0.5; % Location of left base of triangle 2
        end
            
        
        % Linear densities in x for two pops.
        xinit1 = params.grid*params.lengthGx*(mu1+sigma1*sqrt(rand(floor(rho*params.N),1)));
        xinit2 = params.grid*params.lengthGx*(mu2+sigma2*(1-sqrt(1-rand(params.N-floor(rho*params.N),1))));
        yinit = params.grid*params.lengthGy*rand(params.N,1);
        Pos = [[xinit1;xinit2],yinit];
        
        % Impose periodic boundary conditions
        Pos = mod(Pos,params.grid*[params.lengthGx,params.lengthGy]);
        
        % Fixed vectors #2.
        Ang1 = 1*pi/4*ones(floor(params.N/2),1);
        Ang2 = 3*pi/4*ones(params.N-floor(params.N/2),1);
        Ang = [Ang1;Ang2];
        
        BinLocations = [0,pi/2+0.1,pi];
        
    elseif strcmp(IC_type,'WedgeReverse')
            % Reverse Triangle wave setup
        if ~isempty(p.Results.InitialProp)
            rho = p.Results.InitialProp;
        else
            rho = 1/2; %Controls proportion of cells in the two gaussians
        end
        if ~isempty(p.Results.InitialSD)
            sigma1 = p.Results.InitialSD(1); % Base width of triangle 1
            mu1 = 0.5-sigma1; % Location of left base of triangle 1
            sigma2 = p.Results.InitialSD(2); % Base width of triangle 2
            mu2 = 0.5; % Location of left base of triangle 2
        else
            sigma1 = 0.5; % Base width of triangle 1
            mu1 = 0.5-sigma1; % Location of left base of triangle 1
            sigma2 = sigma1; % Base width of triangle 2
            mu2 = 0.5; % Location of left base of triangle 2
        end
        
            % Linear densities in x for two pops.
        xinit1 = mu1+sigma1*(1-sqrt(1-rand(params.N-floor(rho*params.N),1)));
        xinit2 = mu2+sigma2*sqrt(rand(floor(rho*params.N),1));
        yinit = params.grid*params.lengthGy*rand(params.N,1);
        Pos = [[xinit1;xinit2],yinit];
        
        % Impose periodic boundary conditions
        Pos = mod(Pos,params.grid*[params.lengthGx,params.lengthGy]);
        
        % Fixed vectors #2.
        Ang1 = 1*pi/4*ones(floor(params.N/2),1);
        Ang2 = 3*pi/4*ones(params.N-floor(params.N/2),1);
        Ang = [Ang1;Ang2];
        
        BinLocations = [0,pi/2,pi];
        
    elseif strcmp(IC_type,'Wedge2')
            % Squarewave setup
        rho = 1/2; %Controls proportion of cells in the two gaussians
        sigma1 = 0.3875; % Base width of triangle 1
        mu1 = 0.5-sigma1; % Location of left base of triangle 1
        sigma2 = sigma1; % Base width of triangle 2
        mu2 = 0.5; % Location of left base of triangle 2
        
            % Linear densities in x for two pops.
        xinit1 = mu1+sigma1*sqrt(rand(floor(rho*params.N),1));
        xinit2 = mu2+sigma2*(1-sqrt(1-rand(params.N-floor(rho*params.N),1)));
        yinit = params.grid*params.lengthGy*rand(params.N,1);
        Pos = [[xinit1;xinit2],yinit];
        
        % Impose periodic boundary conditions
        Pos = mod(Pos,params.grid*[params.lengthGx,params.lengthGy]);
        
        % Fixed vectors #2.
        Ang1 = 1*pi/4*ones(floor(params.N/2),1);
        Ang2 = 3*pi/4*ones(params.N-floor(params.N/2),1);
        Ang = [Ang1;Ang2];
        
        BinLocations = [0,pi/2,pi];
        
    elseif strcmp(IC_type,'Wedge3')
            % Squarewave setup
        rho = 1/2; %Controls proportion of cells in the two gaussians
        sigma1 = 0.2750; % Base width of triangle 1
        mu1 = 0.5-sigma1; % Location of left base of triangle 1
        sigma2 = sigma1; % Base width of triangle 2
        mu2 = 0.5; % Location of left base of triangle 2
        
            % Linear densities in x for two pops.
        xinit1 = mu1+sigma1*sqrt(rand(floor(rho*params.N),1));
        xinit2 = mu2+sigma2*(1-sqrt(1-rand(params.N-floor(rho*params.N),1)));
        yinit = params.grid*params.lengthGy*rand(params.N,1);
        Pos = [[xinit1;xinit2],yinit];
        
        % Impose periodic boundary conditions
        Pos = mod(Pos,params.grid*[params.lengthGx,params.lengthGy]);
        
        % Fixed vectors #2.
        Ang1 = 1*pi/4*ones(floor(params.N/2),1);
        Ang2 = 3*pi/4*ones(params.N-floor(params.N/2),1);
        Ang = [Ang1;Ang2];
        
        BinLocations = [0,pi/2,pi];
        
    elseif strcmp(IC_type,'Wedge4')
            % Squarewave setup
        rho = 1/2; %Controls proportion of cells in the two gaussians
        sigma1 = 0.1625 ; % Base width of triangle 1
        mu1 = 0.5-sigma1; % Location of left base of triangle 1
        sigma2 = sigma1; % Base width of triangle 2
        mu2 = 0.5; % Location of left base of triangle 2
        
            % Linear densities in x for two pops.
        xinit1 = mu1+sigma1*sqrt(rand(floor(rho*params.N),1));
        xinit2 = mu2+sigma2*(1-sqrt(1-rand(params.N-floor(rho*params.N),1)));
        yinit = params.grid*params.lengthGy*rand(params.N,1);
        Pos = [[xinit1;xinit2],yinit];
        
        % Impose periodic boundary conditions
        Pos = mod(Pos,params.grid*[params.lengthGx,params.lengthGy]);
        
        % Fixed vectors #2.
        Ang1 = 1*pi/4*ones(floor(params.N/2),1);
        Ang2 = 3*pi/4*ones(params.N-floor(params.N/2),1);
        Ang = [Ang1;Ang2];
        
        BinLocations = [0,pi/2,pi];
        
    elseif strcmp(IC_type,'Wedge5')
            % Squarewave setup
        rho = 1/2; %Controls proportion of cells in the two gaussians
        sigma1 = 0.05; % Base width of triangle 1
        mu1 = 0.5-sigma1; % Location of left base of triangle 1
        sigma2 = sigma1; % Base width of triangle 2
        mu2 = 0.5; % Location of left base of triangle 2
        
            % Linear densities in x for two pops.
        xinit1 = mu1+sigma1*sqrt(rand(floor(rho*params.N),1));
        xinit2 = mu2+sigma2*(1-sqrt(1-rand(params.N-floor(rho*params.N),1)));
        yinit = params.grid*params.lengthGy*rand(params.N,1);
        Pos = [[xinit1;xinit2],yinit];
        
        % Impose periodic boundary conditions
        Pos = mod(Pos,params.grid*[params.lengthGx,params.lengthGy]);
        
        % Fixed vectors #2.
        Ang1 = 1*pi/4*ones(floor(params.N/2),1);
        Ang2 = 3*pi/4*ones(params.N-floor(params.N/2),1);
        Ang = [Ang1;Ang2];
        
        BinLocations = [0,pi/2,pi];
        
    end

end

% alignMCwavetest.m
%     % Uniform initial conditions
% %     prop1 = 0.8; %Proportion of cells in main cone
% %     prop2 = 1-prop1; %Proportion of cells in secondary cone
% %     rho1 = 2/3; %Proportion of cells in main cone moving with Ang1
% %     rho2 = 1/3; %Analogue of rho1 for secondary code
% %     xinit = params.grid*params.lengthG*rand(params.N,1);
% %     yinit = params.grid*params.lengthG*rand(params.N,1);
% %     Pos = [xinit,yinit];
% 
%     % Gaussian setup
% %     rho = 1/2; %Controls proportion of cells in the two gaussians
% %     mu1 = [0.25 0.5];
% %     sigma1 = 0.01*eye(2);
% %     mu2 = [0.75 0.5];
% %     sigma2 = 0.01*eye(2);
%     
%     % multivariate
% %     Pos1 = mvnrnd(mu1,sigma1,floor(rho*params.N));
% %     Pos2 = mvnrnd(mu2,sigma2,params.N-floor(rho*params.N));
% %     Pos = [Pos1;Pos2];
%     
%         % univariate in x-direction, uniform in y-direction
% %     Pos1 = mu1(1) + sigma1(1,1)*randn(floor(rho*params.N),1);
% %     Pos2 = mu2(1) + sigma2(1,1)*randn(params.N-floor(rho*params.N),1);
% %     Pos = [Pos1 rand(floor(rho*params.N),1)*params.grid*params.lengthG ...
% %         ;Pos2 rand(params.N-floor(rho*params.N),1)*params.grid*params.lengthG];
%     
%     % Impose periodic boundary conditions
% %     maxI = max(ceil(abs((Pos-0.5)))/(params.grid*params.lengthG),[],'all');
% %     for i = 1:maxI
% %         Pos = Pos + (Pos<0).*params.grid*params.lengthG - (Pos>params.grid*params.lengthG).*params.grid*params.lengthG;
% %     end
%     
%     Pos = mod(Pos,params.grid*params.lengthG);
%     
%     % Plots initial conditions to check periodicity at boundaries.
%         % multivariate
% %     Z1 = mvnpdf(Pos(1:length(Pos1),:),mu1,sigma1);
% %     Z2 = mvnpdf(Pos(length(Pos1)+1:end,:),mu2,sigma2);
%         % univariate in x-direction, uniform in y-direction
% %     Z1 = normpdf(Pos1,mu1(1),sigma1(1,1));
% %     Z2 = normpdf(Pos2,mu2(1),sigma2(1,1));
% 
% %     scatter3(Pos(1:length(Pos1),1),Pos(1:length(Pos1),2),Z1,'r')
% %     hold on
% %     scatter3(Pos(length(Pos1)+1:end,1),Pos(length(Pos1)+1:end,2),Z2,'b')
% %     pause(2)
% %     waitforbuttonpress
%     
% %     xinit1 = params.grid*params.lengthG/2*rand(params.N/2,1);
% %     xinit2 = params.grid*params.lengthG/2*(1+rand(params.N/2,1));
% %     yinit = params.grid*params.lengthG*rand(params.N,1);
% %     Pos = [[xinit1;xinit2],yinit];
%     
%     % Initialize orientation vectors. Set to be uniformly distributed.
% %     Ang = 2*pi*rand(params.N,1)-pi;
% 
%     % Uniform over small angle interval
% %     Ang1 = pi/3*rand(floor(prop1*rho1*params.N),1) - pi/6;
% %     Ang2 = pi/3*rand(ceil(params.N*prop1)-floor(prop1*rho1*params.N),1)+5*pi/6;
% %     Ang3 = 2*pi/3*rand(floor(prop2*rho2*params.N),1)+pi/6;
% %     Ang4 = 2*pi/3*rand(ceil(prop2*params.N)-floor(prop2*rho2*params.N),1)+7*pi/6;
% %     Ang = cat(1,Ang1,Ang2,Ang3,Ang4)-5*pi/6;
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %% alignmentMCParallelogram.m
% 
% 
% 
% 
%       % Simple Uniform;
% %     xinit = params.grid*params.lengthG*rand(params.N,1);
% %     yinit = params.grid*params.lengthG*rand(params.N,1);
% %     Pos = [xinit,yinit];
% 
% 
%     % Wave setup.
% %     rho = 4/7; %Controls proportion of cells in wave
% %     xinit1 = params.grid*params.lengthG/4*rand(floor(rho*params.N),1);
% %     xinit2 = params.grid*params.lengthG*(2/3+1/3*rand(params.N-floor(rho*params.N),1));
% %     yinit = params.grid*params.lengthG*rand(params.N,1);
% %     Pos = [[xinit1;xinit2],yinit];
%     
%         % multivariate
% %     Pos1 = mvnrnd(mu1,sigma1,floor(rho*params.N));
% %     Pos2 = mvnrnd(mu2,sigma2,params.N-floor(rho*params.N));
% %     Pos = [Pos1;Pos2];
    
%     % Impose periodic boundary conditions
% %     maxI = max(ceil(abs((Pos-0.5)))/(params.grid*params.lengthG),[],'all');
% %     for i = 1:maxI
% %         Pos = Pos + (Pos<0).*params.grid*params.lengthG - (Pos>params.grid*params.lengthG).*params.grid*params.lengthG;
% %     end
%         % multivariate
% %     Z1 = mvnpdf(Pos(1:length(Pos1),:),mu1,sigma1);
% %     Z2 = mvnpdf(Pos(length(Pos1)+1:end,:),mu2,sigma2);
%         % univariate in x-direction, uniform in y-direction
% %     Z1 = normpdf(Pos1,mu1(1),sigma1(1,1));
% %     Z2 = normpdf(Pos2,mu2(1),sigma2(1,1));
%     
%     Pos = mod(Pos,params.grid*params.lengthG);
%     
%     % Plots initial conditions to check periodicity at boundaries.
% %     scatter3(Pos(1:length(Pos1),1),Pos(1:length(Pos1),2),Z1,'r')
% %     hold on
% %     scatter3(Pos(length(Pos1)+1:end,1),Pos(length(Pos1)+1:end,2),Z2,'b')
% %     pause(2)
% %     waitforbuttonpress
%     
% %     xinit1 = params.grid*params.lengthG/2*rand(params.N/2,1);
% %     xinit2 = params.grid*params.lengthG/2*(1+rand(params.N/2,1));
% %     yinit = params.grid*params.lengthG*rand(params.N,1);
% %     Pos = [[xinit1;xinit2],yinit];
%     
%     % Initialize orientation vectors. Set to be uniformly distributed.
% %     Ang = 2*pi*rand(params.N,1)-pi;
% 
%     % Uniform over small angle interval
% %     Ang1 = pi/3*rand(params.N/2,1);
% %     Ang2 = pi/3*rand(params.N/2,1)+2*pi/3;
% %     Ang = [Ang1;Ang2];
% 

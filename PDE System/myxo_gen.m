function [] = myxo_gen(Num,kappa,T,results_folder)

% clear all
% clear global


% we're solving the system
%
% d_t p1 + c1 d_x p1(1 + c2 p2) = 0
% d_t p2 - c1 d_x p2(1 + c2 p1) = 0
%

% global ilim
% ilim = 0; % 0 = low-order, 1 = limited second order

nlev = 3; % mesh refinement level


tmesout = 10;   % time-interval for messaging
tmes = 0;


LL = 400;   % domain

% parameters in the model
zeta1 = [cos(pi/4), sin(pi/4)];
zeta2 = [cos(3*pi/4), sin(3*pi/4)];
c1 = 1; %zeta1(1)*vbar;
c2 = (kappa/Num)*LL^2; %abs([-zeta1(2), zeta1(1)] * zeta2')*length^2/2;
% Num = 2856;%2856*16;45700

% K = additional diffusion in the definition of the flux
K = c1*(1 + c2);

Delt = 0.1;   % output step

% time-step of simulation
% number of space points
dt = 0.1e-2 / 2^nlev;
N = 400 * 2^nlev * 10;

% time of simulation
% T should be an integer
% T = 50; 
ic1 = 4; %initial conditions: (only 4 and 5 set up for general N)
         %   1 - square (4000 cells per band)
         %   2 - normal (2000 cells per band)
         %   3 - smoothed square (4000 cells per band)
         %   4 - normal+background density (1000 cells per band)
         %   5 - one normal band+background density (? cells per orientation)
         
tspan = linspace(0, T, T/Delt+1);


dx = LL / N;
x = zeros(N, 1);
for k=1:N
    x(k) = (k-1)*dx;
end

alpha = K*dt / dx
addl_diff = K/2*dx


% initial conditions
y0 = zeros(2*N,1);
% IC

if ic1 == 1
    y0(1:N) = y0(1:N) + 5*(40/400) * (x >= 100) .* (x <= 200);
    y0((N+1):2*N) = y0((N+1):2*N) + 5*(40/400) * (x >= 200) .* (x <= 300);
end
if ic1 == 2
    y0(1:N) = y0(1:N) +  (500/400).*normpdf(x,150,25); %exp(-(x - 150).^2 /625 );
    y0((N+1):2*N) = y0((N+1):2*N) + (500/400).*normpdf(x,250,25); %exp(-(x - 250).^2 /625);
    
end
if ic1 == 3
    fn1 = @(x) exp(-((x - 150)/50).^8);
    fn2 = @(x) exp(-((x - 250)/50).^8);
    n1 = integral(fn1,0,400);
    n2 = integral(fn2,0,400);
    y0(1:N) = y0(1:N) +  (20000/(400*n1)).*exp(-((x - 150)/50).^8  );
    y0((N+1):2*N) = y0((N+1):2*N) + (20000/(400*n2)).*exp(-((x - 250)/50).^8 );
end
if ic1 == 4
    y0(1:N) = y0(1:N) +  1/8*(Num/2).*normpdf(x,150,25)+7/8*(Num/2/400); %0.25*exp(-(x - 150).^2 /(2*625) ) + 0.75;
    y0((N+1):2*N) = y0((N+1):2*N) + 1/8*(Num/2).*normpdf(x,250,25)+7/8*(Num/2/400); %0.25*exp(-(x - 250).^2 /(2*625))+0.75;
    y0 = y0/LL;% creates density of a slice in the x-direction (L_y = 400)
end
if ic1 == 5
    y0(1:N) = y0(1:N) +  1/15*(Num).*normpdf(x,150,25)+7/15*(Num/LL); %0.25*exp(-(x - 150).^2 /(2*625) ) + 0.75;
    y0((N+1):2*N) = y0((N+1):2*N) + 7/15*(Num/LL); %0.25*exp(-(x - 250).^2 /(2*625))+0.75;
    y0 = y0/LL; % creates density of a slice in the x-direction (L_y = 400)
end

pltic = 0;
if pltic > 0
    figure(3);
    clf;
    plot(x,y0(1:N),'b-', x,y0((N+1):2*N),'r-');
    axis([0 LL 0 1]);
    title('IC');
    pause(0.2);
end

% Integrator in time
% use_eul = 1 to use Euler
use_eul = 1;

if use_eul > 0
    numsteps = Delt/dt;
    indx = 1;
    tt(indx) = 0;
    y(indx,:) = y0;
    soln = y0;
    time = 0;

    for kk=0:Delt:T
        for jj=1:numsteps
          time = time + dt; 
          [yprime, tmes] = rhsm(time, soln, dx,c1,c2,N,K,tmes,tmesout);
          
          % second-order time-stepping
          soln1 = soln + dt*yprime;
          [yprime, tmes] = rhsm(time, soln1, dx,c1,c2,N,K,tmes,tmesout);
          soln = (soln + soln1 + dt*yprime)/2;
          
          % Euler time-stepping
          %soln = soln + dt*yprime;
        end

        indx = indx+1;
        y(indx,:) = soln;
        tt(indx)   = time;
    end
else
    rhsm_handle = @(t,y)rhsm(t, y, dx,c1,c2,N,K,tmes,tmesout);
    [tt,y] = ode45(@rhsm_handle, tspan, y0);
end


save([results_folder '/mean field data kappa = ' num2str(kappa) ', N = ' num2str(Num) ' alpha = ' num2str(alpha) ' addl_diff = ' num2str(addl_diff) '.mat'],'Num','kappa','x','tt','y')


%% plotting

v = VideoWriter([results_folder '/mean field data kappa = ' num2str(kappa) ', N = ' num2str(Num) ' alpha = ' num2str(alpha) ' addl_diff = ' num2str(addl_diff)]);
v.FrameRate = 20;
open(v)

sizey = size(y);

clear h
h = figure();
clf;


for kk=1:sizey(1)
        
    plot(x,y(kk,1:N),'b-', x,y(kk,(N+1):end),'r-','Linewidth', 2);
    ylim([0 max(y(:))+0.2*max(y(:))]);
    time = tt(kk);
    title(strcat("time=",num2str(time)));
    
    M = getframe(gcf);
    writeVideo(v,M)
    xlabel('position, x')
    ylabel('linear cell density')
    
end

hold on
plot(x,y(1,1:N),'k-', x,y(1,(N+1):end),'k-','Linewidth', 2);

M = getframe(gcf);
writeVideo(v,M)
close(v)

title(['T = ' num2str(T) ' with N = ' num2str(Num) ', length = ' num2str((2*c2)^0.5) ', kappa = ' num2str(Num*c2/LL^2)])
legend('right-moving wave','left-moving wave','right-moving initial profile',...
        'left-moving initial profile','Location','south')
    
savefig([results_folder '/Final and initial profiles kappa = ' num2str(kappa) ' N = ' num2str(Num) ' alpha = ' num2str(alpha) ' addl_diff = ' num2str(addl_diff)  '.fig'])
saveas(h,[results_folder '/Final and initial profiles kappa = ' num2str(kappa) ' N = ' num2str(Num) ' alpha = ' num2str(alpha) ' addl_diff = ' num2str(addl_diff) '.png'],'png')




end


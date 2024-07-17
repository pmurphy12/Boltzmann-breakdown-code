
clear all

load run.mat




sizey = size(y);

figure()
clf;


% y=400*y; % normalization to linear number density    
for kk=1:sizey(1)
        
    plot(x,y(kk,1:N),'b-', x,y(kk,(N+1):end),'r-','Linewidth', 2);
    ylim([0 max(y(:))+0.2*max(y(:))]);
    time = tt(kk);
    title(strcat("time=",num2str(time)));
    pause(0.0001);
    
end

hold on
plot(x,y(1,1:N),'k-', x,y(1,(N+1):end),'k-','Linewidth', 2);
xlabel('position, x')
ylabel('linear cell density')
title(['T=50 with N = ' num2str(Num) ', length = ' num2str((2*c2)^0.5) ', kappa = ' num2str(Num*c2/LL^2)])
% title(['T=0.125 with p0 = ' num2str(Num) ', kappa = ' num2str(c2)])



% Wave Speed Calc

SampleRate = 100;
L = length(data);
xPosMean = zeros(L,1);
ang = 3*pi/4;
EPS = 10^-8;
frac = 0.2;

for i = 1:L
    pos = data{i}.pos;
    I = data{i}.o == ang;
%     I = data{i}.o <= ang+EPS & data{i}.o >= ang-EPS;
    
    xPos = sort(pos(I,1));
    
    xPos = xPos(1:floor(length(xPos)*0.2));
    
    xPosMean(i) = mean(xPos);
    
end

figure;
plot(linspace(0,params.tend,101),xPosMean)
figure;
plot(linspace(0,params.tend,100),diff(xPosMean)/(params.dt*SampleRate))



if not(exist("dataCat") == 1)

    load('Concatenated rod info.mat', 'dataCat')

end
nEdges = 401;
I = dataCat{end}.o == pi/4;
x = dataCat{end}.pos(I,1);
y = dataCat{end}.pos(I,2);
bins = linspace(0,400,nEdges);
KS = zeros(length(bins)-1,1);
P = KS;

% unifSample = 400*rand(length(y),1);

for i = 1:(length(bins)-1)

    Itemp = (x >= bins(i)) & (x <= bins(i+1));
    S = y(Itemp);
    unifSample = 400*rand(length(S),1);
    [KS(i),P(i)] = kstest2(S,unifSample,'Alpha',0.05);

end

reject_percentage = sum((KS == 1))/length(KS)

figure;
b = bins(1:end-1)/400;
scatter(b(KS == 0),KS(KS == 0),'o')
hold on
scatter(b(KS == 1),KS(KS == 1),'o','r')
hold off
xlabel('x')
ylabel('Keep null hypothesis (0) or reject (1)')
title('')



% figure;
% surf(X,Y,0.5*D/(sum(D,'all')*(Y(2,1)-Y(1,1))*(X(1,2)-X(1,1))),'EdgeColor','none')


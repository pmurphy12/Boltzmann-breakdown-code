
if not(exist("dataCat") == 1)

    load('Concatenated rod info.mat', 'dataCat')

end
I = dataCat{end}.o == pi/4;
x = dataCat{end}.pos(I,1);
y = dataCat{end}.pos(I,2);
[bdw,D,X,Y] = kde2d([x,y],1048,[0 0], [400 400]);

KLD = zeros(size(D,2),1);
dy = Y(2,1)-Y(1,1);

for i = 1:size(D,2)

    P = D(:,i)/(sum(D(:,i))*dy);
    Q = ones(size(D,2),1)/(size(D,2)*dy);
    KLD(i) = KLtest_continuous(P,Q,dy);

end

figure;
plot(Y(:,1),KLD)
xlabel('x')
ylabel('KL value')

figure;
surf(X,Y,0.5*D/(sum(D,'all')*(Y(2,1)-Y(1,1))*(X(1,2)-X(1,1))),'EdgeColor','none')


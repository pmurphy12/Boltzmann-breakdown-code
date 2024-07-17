

Ns = [100:100:10000];
simN = 50;
meanMinD = zeros(size(Ns));

for i = 1:length(Ns)
    
    temp = cell(1,simN);
    
    for ii = 1:simN
        x = 0.25*rand(Ns(i),1);
        y = rand(Ns(i),1);
        D = pdist2([x(1:Ns(i)/2),y(1:Ns(i)/2)],[x(Ns(i)/2+1:end),y(Ns(i)/2+1:end)],'euclidean','Smallest',2);
%         temp(ii) = mean(D(2,:));
        temp{ii} = D(2,:);
    end
    
    meanMinD(i) = mean(cat(2,temp{:}));
    
end

figure; plot(Ns,meanMinD)
axis([min(Ns) max(Ns) 0 max(meanMinD)])

g = fittype('a*x+b');
[curve, goodness] = fit( log(Ns)', log(meanMinD'), g);
curve

hold on
plot(Ns,exp(curve(log(Ns))))

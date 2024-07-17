p11 = 0.5;
p12 = 0.5;
p21 = 0.3;
p22 = 0.7;

k = 1000;
c1 = zeros(k,2);
c2 = zeros(k,2);
for i = 1:k
n1 = randi([0 3],1);
n2 = randi([2 7],1);
c1(i,1) = binornd(n1,p11);
c1(i,2) = n1 - c1(i,1);
c2(i,1) = binornd(n2,p21);
c2(i,2) = n2 - c2(i,1);
end

% mean(c1(:,1).*c1(:,2)./(sum(c1,2))./(sum(c1,2)-1),'omitnan')-mean(c1(:,1)./(sum(c1,2)),'omitnan')*mean((c1(:,2))./(sum(c1,2)),'omitnan')
% mean(c2(:,1).*c2(:,2)./(sum(c2,2))./(sum(c2,2)-1),'omitnan')-mean(c2(:,1)./(sum(c2,2)),'omitnan')*mean((c2(:,2))./(sum(c2,2)),'omitnan')

(mean(c1(:,1).*c1(:,2),'omitnan')-mean(c1(:,1),'omitnan')*mean((c1(:,2)),'omitnan'))./(mean(c1(:,1),'omitnan')*mean((c1(:,2)),'omitnan'))
(mean(c2(:,1).*c2(:,2),'omitnan')-mean(c2(:,1),'omitnan')*mean((c2(:,2)),'omitnan'))./(mean(c2(:,1),'omitnan')*mean((c2(:,2)),'omitnan'))





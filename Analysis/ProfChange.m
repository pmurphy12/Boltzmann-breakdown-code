function [Stats] = ProfChange(data,params,varargin)
%   Calculate interquartile range, skewness, and change in front and back
%   of support above EPS for profiles.

    PLOT = false;

    ANG = unique(data{1}.o);
    Stats.EPS = 10^-3;
    
    IQR = zeros(length(ANG),length(data));
    Skew = zeros(length(ANG),length(data));
    SupLeft = zeros(length(ANG),length(data));
    SupRight = zeros(length(ANG),length(data));
    m = zeros(length(ANG),length(data));
    
    for i = 1:length(data)
        k = 1;
        for j = ANG'
            I = (data{i}.o <= j+Stats.EPS) & (data{i}.o >= j-Stats.EPS);
            
            IQR(k,i) = iqr(data{i}.pos(I,1));
            Skew(k,i) = skewness(data{i}.pos(I,1));
            
            h = histcounts(data{i}.pos(I,1),0:params.lengthGx/500:params.lengthGx);
            yL = cumsum(h);
            SupLeft(k,i) = find(yL/sum(I)>Stats.EPS,1,'first');
            SupRight(k,i) = find(yL/sum(I)>1-Stats.EPS,1,'first');
            
            m(k,i) = mean(data{i}.pos(I,1));
            
            k = k+1;
        end
    end
    
    l = 1;
    for i = [1 length(data)]
        k = 1;
        for j = ANG'
            I = data{i}.o == j;
            dataS{k,l} = data{i}.pos(I,1);
            Mode(k,l) = findMode(data{i}.pos(I,1));
            Median(k,l) = median(data{i}.pos(I,1));
            
            k = k+1;
        end
        l = l+1;
    end
    
    
    
    skQ = @(Q) (Q(3)+Q(1)-2*Q(2))/(Q(3)-Q(1)); % Bowley's method of skewness
    skp1 = @(x) (mean(x) - findMode(x))/std(x); % Pearson's first skewness coeff
    skp2 = @(x) 3*(mean(x) - median(x))/std(x); % Pearson's second skewness coeff
    abCalc1 = @(S,m) abs(S(:,1) - m(:,1));
    abCalcEnd = @(S,m) abs(S(:,end) - m(:,end));
    Rel = @(x1,xend) (xend-x1)./x1;
    stdL = @(x,m) sqrt(1./(length(x(x<m))-1).*sum((x(x<m)-m).^2));
    stdR = @(x,m) sqrt(1./(length(x(x>m))-1).*sum((x(x>m)-m).^2));
    
    abRatio1a = @(f,SL,SR,m) f(SL,m)./f(SR,m);
    abRatio2a = @(f,SL,SR,m) f(SL,m)./(f(SL,m)+f(SR,m));
    
    Q1 = quantile(dataS{1,1},[0.25,0.5,0.75]);
    Qend = quantile(dataS{1,2},[0.25,0.5,0.75]);
    
    Stats.Bowley1 = skQ(Q1);
    Stats.BowleyEnd = skQ(Qend);
    Stats.BowleyChange = Stats.BowleyEnd - Stats.Bowley1;
    Stats.BowleyRelChange = Rel(Stats.Bowley1,Stats.BowleyEnd);
    
    Stats.Pearson11 = skp1(dataS{1,1});
    Stats.Pearson1end = skp1(dataS{1,2});
    Stats.Pearson1Change = Stats.Pearson1end - Stats.Pearson11;
    Stats.Pearson1RelChange = Rel(Stats.Pearson11,Stats.Pearson1end);
    
    Stats.Pearson21 = skp2(dataS{1,1});
    Stats.Pearson2end = skp2(dataS{1,2});
    Stats.Pearson2Change = Stats.Pearson2end - Stats.Pearson21;
    Stats.Pearson2RelChange = Rel(Stats.Pearson21,Stats.Pearson2end);
    
    Stats.abRatio1a1 = abRatio1a(abCalc1,SupLeft,SupRight,m);
    Stats.abRatio1aEnd = abRatio1a(abCalcEnd,SupLeft,SupRight,m);
    Stats.abRatio1aChange = Stats.abRatio1aEnd - Stats.abRatio1a1;
    Stats.abRatio1aRelChange = Rel(Stats.abRatio1a1,Stats.abRatio1aEnd);
    
    Stats.abRatio1b1 = 1./Stats.abRatio1a1;
    Stats.abRatio1bEnd = 1./Stats.abRatio1aEnd;
    Stats.abRatio1bChange = Stats.abRatio1bEnd - Stats.abRatio1b1;
    Stats.abRatio1bRelChange = Rel(Stats.abRatio1b1,Stats.abRatio1bEnd);
    
    Stats.abRatio2a1 = abRatio2a(abCalc1,SupLeft,SupRight,m);
    Stats.abRatio2aEnd = abRatio2a(abCalcEnd,SupLeft,SupRight,m);
    Stats.abRatio2aChange = Stats.abRatio2aEnd - Stats.abRatio2a1;
    Stats.abRatio2aRelChange = Rel(Stats.abRatio2a1,Stats.abRatio2aEnd);
    
    Stats.abRatio2b1 = 1-Stats.abRatio2a1;
    Stats.abRatio2bEnd = 1-Stats.abRatio2aEnd;
    Stats.abRatio2bChange = -Stats.abRatio2aChange;
    Stats.abRatio2bRelChange = Rel(Stats.abRatio2b1,Stats.abRatio2bEnd);
    
    
    Stats.SupRChange = (SupRight(:,end) - m(:,end)) - (SupRight(:,1) - m(:,1));
    Stats.SupLChange = (SupLeft(:,end) - m(:,end)) - (SupLeft(:,1) - m(:,1));

    Stats.LeftSD1 = stdL(dataS{1,1},m(1,1));
    Stats.LeftSDend = stdL(dataS{1,2},m(1,end));
    Stats.LeftSDChange = Stats.LeftSDend - Stats.LeftSD1;
    Stats.LeftSDRelChange = Rel(Stats.LeftSD1,Stats.LeftSDend);
    
    Stats.RightSD1 = stdR(dataS{1,1},m(1,1));
    Stats.RightSDend = stdR(dataS{1,2},m(1,end));
    Stats.RightSDChange = Stats.RightSDend - Stats.RightSD1;
    Stats.RightSDRelChange = Rel(Stats.RightSD1,Stats.RightSDend);

    Stats.ANG = ANG;
    
    Stats.Skew = Skew;
    Stats.SkewChange = Skew(:,end) - Skew(:,1);
    Stats.IQR = IQR;
    Stats.IQRchange = IQR(:,end) - IQR(:,1);
    Stats.SupLeft = SupLeft;
    Stats.SupRight = SupRight;
    
    
%     F = fieldnames(StatsSum);
% 
%     for i = 1:length(F)-1
%     A.(F{i}) = mean(cat(2,StatsSum(:).(F{i})),2);
%     end
%     s
    
    if PLOT
        for i = 1:length(ANG)
            figure; plot(Skew(i,:))
            title(['Angle ' num2str(ANG(i))])
            xlabel('time')
            ylabel('Skew')

            figure; plot(IQR(i,:))
            title(['Angle ' num2str(ANG(i))])
            xlabel('time')
            ylabel('IQR')

            figure; plot(SupLeft(i,:))
            title(['Angle ' num2str(ANG(i))])
            xlabel('time')
            ylabel('SupLeft')

            figure; plot(SupRight(i,:))
            title(['Angle ' num2str(ANG(i))])
            xlabel('time')
            ylabel('SupRight')
        end
    end

end


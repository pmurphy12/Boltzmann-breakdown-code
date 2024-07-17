function [statOut] = statAnalysis(data)
%STATANALYSIS Summary of this function goes here
%   Detailed explanation goes here
    
    L = length(data);
    cMean = zeros(L,2);
    cVar = zeros(L,2);
    cSTD = zeros(L,2);
    cDisp = zeros(L,2);
    Support = zeros(L,2);
    Mean = zeros(L,2);
    Var = zeros(L,2);
    Skew = zeros(L,2);
    

    parfor i = 1:L
        Pos = data{i}.pos;
        Ang = data{i}.o;
        ang = unique(Ang);
        pos1 = Pos(Ang == ang(1),1);
        pos2 = Pos(Ang == ang(2),1);
        pos = {pos1, pos2};
        
        for j = 1:2
            [cMean(i,j),cVar(i,j),cSTD(i,j),cDisp(i,j)] = circStatsSample(pos{j},1);
            
            %Shift so mean is centered at 0.5
            pos{j} = mod(pos{j} - cMean(i,j)+0.5,1);
        
            Support(i,j) = max(pos{j})-min(pos{j});
        
            Mean(i,j) = mean(pos{j});
        
            Var(i,j) = var(pos{j});
        
            sigma = sqrt(Var(i,j));
        
            Skew(i,j) = mean(((pos{j}-Mean(i,j))./sigma).^3);
        end


    end
    
    statOut = struct('cMean',cMean,'cVar',cVar,'cSTD',cSTD,'cDisp',cDisp,...
                        'Support',Support,'Mean',Mean,'Var',Var,'Skew',Skew);
    
end


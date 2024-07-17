function [] = statsComp(Stats1,Stats2,Stats3,Name)


    F = fieldnames(Stats1);
%     F = intersect(F1,F2);
%     A.('SupLChange') = A.('SupLChange')/500;
    F = setdiff(F,{'ANG','EPS','SupLeft','SupRight','Type','IQR','Skew'});
%     c = categorical({'Agent Simulation','Numerical Solution'});

    for i = 1:length(F)
        if size(cat(2,Stats1.(F{i})),1) == 2
            a = cat(2,Stats1(:).(F{i}));
            a = a(1,:);
            b = cat(2,Stats2(:).(F{i}));
            b = b(1,:);
            c = cat(2,Stats3(:).(F{i}));
            c = c(1,:);
        else
            a = cat(2,Stats1(:).(F{i}));
            b = cat(2,Stats2(:).(F{i}));
            c = cat(2,Stats3(:).(F{i}));
        end
        
        h = figure;
        scatter(0.5*ones(length(a),1),a)
        hold on
        scatter(0.05*ones(length(b),1),b)
        scatter(0.005*ones(length(c),1),c)
        title(F{i})
        saveas(h,[Name ' ' F{i}],'png')
        
        h1 = figure;
        loglog(0.5*ones(length(a),1),abs(a),'o')
        hold on
        loglog(0.05*ones(length(b),1),abs(b),'o')
        loglog(0.005*ones(length(c),1),abs(c),'o')
        title(F{i})
        saveas(h1,[Name ' log scale ' F{i}],'png')
        
    end



end
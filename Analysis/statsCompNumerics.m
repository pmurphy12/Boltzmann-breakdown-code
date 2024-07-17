function [] = statsCompNumerics(Stats,Name)


    F = fieldnames(Stats);
    F = F(contains(F,'Change'));
    
    kappa = cat(2,Stats(:).('kappa'));
    
    for i = 1:length(F)
        a = cat(2,Stats(:).(F{i}));
        
        h = figure;
        scatter(kappa,a)
        title(F{i})
        savefig(h,[Name ' ' F{i}])
        saveas(h,[Name ' ' F{i}],'png')
        
        h1 = figure;
        loglog(kappa,abs(a),'o')
        title(F{i})
        savefig(h1,[Name ' ' F{i}])
        saveas(h1,[Name ' log scale ' F{i}],'png')
        
    end



end
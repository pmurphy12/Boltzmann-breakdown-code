function [] = plotFeatureKappa(Tcomp,feature,time,yshift)

% Tcomp is table of feature comparisons between mean field and agent based model
% feature is specific feature to compare against.
% ffeatures: 'peakLocation', 'peakAmplitude','skewness'

    SAVEFIG = true;
    save_folder = './';
    
    f = unique(Tcomp.feature);
    if ~any(strcmp(f,feature))
        error('unused feature. check spelling or add desired feature to previous code profileAnalysisRun.m')
    end
    
    T = Tcomp(strcmp(Tcomp.feature,feature),:);
    T = T(T.time == time,:);
    T = T(T.yshift == yshift,:);
    
    h = figure;
    hold on
    
    plot(T.('agent kappa'),T.('mean-field kappa'))
    ylabel('mean-field kappa')
    xlabel('agent sim kappa')
    title(['comparison for feature ' feature ' at time ' num2str(time) ' with y-shift ' num2str(yshift)])
    
    if SAVEFIG
        savefig([save_folder 'abs-mf kappa comparison for feature ' feature ' at time ' num2str(time) ' with y-shift ' num2str(yshift)])
        saveas(h,[save_folder 'abs-mf kappa comparison for feature ' feature ' at time ' num2str(time) ' with y-shift ' num2str(yshift) '.png'],'png')
    end
    
    
end
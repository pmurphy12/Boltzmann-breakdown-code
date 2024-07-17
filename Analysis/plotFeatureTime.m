function [] = plotFeatureTime(T_abs,T_mf,feature,kappa,yshift)

% Tcomp is table of feature comparisons between mean field and agent based model
% feature is specific feature to compare against.
% ffeatures: 'peakLocation', 'peakAmplitude','skewness'

    SAVEFIG = true;
    save_folder = './';
    
%     f = unique(T_abs.feature);
%     if ~any(strcmp(f,feature))
%         error('unused feature. check spelling or add desired feature to previous code profileAnalysisRun.m')
%     end
    
    Ta = T_abs;
    Ta.kappa = round(Ta.kappa*1000)/1000;
    Ta = Ta(Ta.kappa == kappa,:);
    Ta = Ta(Ta.yshift == yshift,:);
    
    Tm = T_mf;
    if sum(Tm.kappa == kappa) == 0
        [~,I] = min(abs(Tm.kappa - kappa));
        kappa2 = Tm.kappa(I)
        Tm = Tm(Tm.kappa == kappa2,:);
    else
        Tm = Tm(Tm.kappa == kappa,:);
    end
    
    if isempty(Ta) || isempty(Tm)
        error('either mean-field data or abs data table is empty based on chosen parameters')
    end
    
    h = figure;
    hold on
    
    plot(Ta.time,Ta.(feature))
    plot(Tm.time,Tm.(feature))
    ylabel(feature)
    xlabel('time')
    legend('agent-based model','mean-field model','Location','northwest')
    
    if SAVEFIG
        savefig([save_folder 'abs and mf time course for feature ' feature ' with kappa ' sprintf('%1.1e',kappa) ' with y-shift ' num2str(yshift) '.fig'])
        saveas(h,[save_folder 'abs and mf time course for feature ' feature ' with kappa ' sprintf('%1.1e',kappa) ' with y-shift ' num2str(yshift) '.png'],'png')
    end
    
    
end
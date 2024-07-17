function [Tcomp] = profileAnalysisComp(T_mf,T_abs)
%profileAnalysisComp.m finds closest mean-field match to agent profiel
%based on given feature(s)
    
    Title = 'full sim';
    SAVEDATA = true;

    T_abs = T_abs(~isnan(T_abs.yshift),:); % removes undefined y-shifts
    
    feats = {'peakLocation','peakAmplitude',...
        'skewness'};
    varNames = {'feature','feature value','time','absolute error',...
        'relative error','agent kappa','mean-field kappa','yshift'};
    
    feature = cell(height(T_abs),length(feats));
    feature_val = zeros(height(T_abs),length(feats));
    time = zeros(height(T_abs),length(feats));
    comp_error = zeros(height(T_abs),length(feats));
    comp_relerror = zeros(height(T_abs),length(feats));
    kappa_abs = zeros(height(T_abs),length(feats));
    kappa_mf = zeros(height(T_abs),length(feats));
    yshift = zeros(height(T_abs),length(feats));
    
    for k = 1:length(feats)
        for i = 1:height(T_abs)
            feature{i,k} = feats{k};
            kappa_abs(i,k) = T_abs.kappa(i);
            yshift(i,k) = T_abs.yshift(i);
            
            time(i,k) = T_abs.time(i);
            [~,I_t_mf] = min(abs(time(i,k)-T_mf.time));
            t_mf = T_mf.time(I_t_mf);
            Ti = T_mf(T_mf.time == t_mf,:);
            
            feature_val(i,k) = T_abs.(feats{k})(i);
            [~,I] = min(abs(feature_val(i,k)-Ti.(feats{k})));
            feat_mf = Ti.(feats{k})(I);
            
            comp_error(i,k) = feature_val(i,k)-feat_mf;
            comp_relerror(i,k) = (feature_val(i,k)-feat_mf)/feature_val(i,k);
            
            kappa_mf(i,k) = Ti.kappa(I);
            
        end
    end
    
    Tcomp = table(feature(:),...
                feature_val(:),...
                time(:),...
                comp_error(:),...
                comp_relerror(:),...
                kappa_abs(:),...
                kappa_mf(:),...
                yshift(:),...
                'VariableNames',varNames);
            
    if SAVEDATA
        save(['Agent profile features comparison to mean field ' Title],'Tcomp','T_abs','T_mf')
    end

end


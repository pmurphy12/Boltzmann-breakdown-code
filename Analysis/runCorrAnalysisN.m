function [time_stats] = runCorrAnalysisN(tSteps,identificationString,varargin)

% Runs correlationAnalysisN.m for the inputed times and saves both spatial
% correlation data for each as well as a spatial average and plots the
% resulting time-course data.
    
    SAVEFIG = true;
    Fact = 1;

    %preallocate
    time_stats.Stats = cell(length(tSteps),1);
    
    time_stats.meanCorr11 = zeros(length(tSteps),1);
    time_stats.meanCorr12 = zeros(length(tSteps),1);
    time_stats.meanCorr22 = zeros(length(tSteps),1);
    time_stats.meanCorrRatio11 = zeros(length(tSteps),1);
    time_stats.meanCorrRatio12 = zeros(length(tSteps),1);
    time_stats.meanCorrRatio22 = zeros(length(tSteps),1);
    time_stats.meanCorr = zeros(length(tSteps),1);
    
    N = 2^6;
    
    %run loop
    for i = 1:length(tSteps)
        [Stats,~,~,~,params] = correlationAnalysisN(tSteps(i),identificationString,'N',N,'Fact',Fact,'SAVEFIG',true);
        
        time_stats.meanCorr11(i) = mean(abs(Stats.independence11(~isinf(Stats.independence11))),'all','omitnan');
        time_stats.meanCorr12(i) = mean(abs(Stats.independence12(~isinf(Stats.independence12))),'all','omitnan');
        time_stats.meanCorr22(i) = mean(Stats.independence22(~isinf(Stats.independence22)),'all','omitnan');
        time_stats.meanCorrRatio11(i) = mean(abs(Stats.independenceRatio11(~isinf(Stats.independenceRatio11))),'all','omitnan');
        time_stats.meanCorrRatio12(i) = mean(abs(Stats.independenceRatio12(~isinf(Stats.independenceRatio12))),'all','omitnan');
        time_stats.meanCorrRatio22(i) = mean(abs(Stats.independenceRatio22(~isinf(Stats.independenceRatio22))),'all','omitnan');
        time_stats.meanCorr(i) = mean(abs(Stats.fixedSpaceCorr(~isinf(Stats.fixedSpaceCorr))),'all','omitnan');
        
        time_stats.Stats{i} = Stats;
    end
    
    save(['Statistical independence test with discretization ' num2str(N)],'time_stats')
    
    h1 = figure;
    hold on;
    t = tSteps*params.dt;
    plot(t,time_stats.meanCorr11)
    plot(t,time_stats.meanCorr12)
    plot(t,time_stats.meanCorr22)
    
    legend('angle 1: mean # of pairs - mean #^2','angles 1&2: mean # of pairs - mean # 1 * mean # 2','angle 2: mean # of pairs - mean #^2','Location','northwest')
    if SAVEFIG
            savefig(['average from local independence test with discretization ' num2str(N)])
            saveas(h1,['average from local independence test with discretization ' num2str(N) '.png'],'png')
    end
    
    h2 = figure;
    hold on;
    plot(t,time_stats.meanCorrRatio11)
    plot(t,time_stats.meanCorrRatio12)
    plot(t,time_stats.meanCorrRatio22)
    axis([min(t) max(t) -1/2 1/2]);
    legend('angle 1: (mean # of pairs - mean #^2)/mean # of pairs','angles 1&2: (mean # of pairs - mean # 1 * mean # 2)/mean # of pairs','angle 2: (mean # of pairs - mean #^2)/mean # of pairs','Location','northwest')
    if SAVEFIG
            savefig(['average from local independence ratio test with discretization ' num2str(N)])
            saveas(h2,['average from local independence ratio test with discretization ' num2str(N) '.png'],'png')
    end
    
    h3 = figure;
    hold on;
    plot(t,time_stats.meanCorr)
    legend('correlation between the 2 angles','Location','northwest')
    if SAVEFIG
            savefig(['average from correlation test with discretization ' num2str(N)])
            saveas(h3,['average from correlation test with discretization ' num2str(N) '.png'],'png')
    end
    

end
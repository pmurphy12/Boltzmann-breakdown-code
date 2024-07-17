function [T_mf] = profileAnalysisRun(identificationString)
% Calculate table of mean-field profile properties. Meant to be run from
% the "Mean Field" folder.
% identificationString = 'kappa'

    Title = 'full sim';
    SAVEFIG = true;
    SAVEDATA = true;
    
    addpath('./mean-field data')

    Files = dir(fullfile('./mean-field data/',['*' identificationString '*.mat']));
    
    varNames = {'time','peakLocation','peakAmplitude','skewness','kappa','Num','wave'};
    
    % extract information from files
    listOfVariables = who('-file',['./mean-field data/' Files(1).name]);
    load(Files(1).name,'tt');
%     tt = 3200;

    time = zeros(length(Files),length(tt));
    kappa = zeros(length(Files),length(tt));
    Num = zeros(length(Files),length(tt));
    peakLocation = zeros(length(Files),length(tt));
    peakAmplitude = zeros(length(Files),length(tt));
    skewness = zeros(length(Files),length(tt));
    wave_direction = cell(length(Files),length(tt));
    
    
    for ii = 1:length(Files)
            S = load(Files(ii).name,listOfVariables{:});
            y = S.y;
            
        parfor jj = 1:length(S.tt)
            time(ii,jj) = S.tt(jj);
            kappa(ii,jj) = S.kappa;
            Num(ii,jj) = S.Num;
            % get statistics for profiles
            [peakLocation(ii,jj),peakAmplitude(ii,jj),skewness(ii,jj)] = profileAnalysis(S.x,y(jj,1:length(S.x)));
            wave_direction{ii,jj} = 'right-moving';
    
        end
    end
    
    
    % Table generation
    
    
    T_mf = table(time(:),...
                peakLocation(:),...
                peakAmplitude(:),...
                skewness(:),...
                kappa(:),...
                Num(:),...
                wave_direction(:),...
                'VariableNames',varNames);
    
    if SAVEDATA
        save(['Mean Field profile features ' Title],'T_mf')
    end
    
%     %plotting in time
%     clear h
%     h = figure;
%     tl = tiledlayout(1,1);
%     ax1 = axes(tl);
%     plot(ax1,t(I),Stats.mean_in_time)
%     ax1.YLim = [-0.2 0.2];
%     xlabel('time')
%     ylabel('mean asymmetry (time)')
%     
%     ax2 = axes(tl);
%     plot(ax2,x(1:end-1),Stats.mean_in_space,'Color',[0.8 0.2 0.2])
%     ax2.XAxisLocation = 'top';
%     ax2.YAxisLocation = 'right';
% %     ax2.YLim = [-0.2 0.2];
%     xlabel('x')
%     ylabel('mean asymmetry (space)')
%     
%     ax2.Color = 'none';
%     ax2.XColor = [0.8 0.2 0.2];
%     ax2.YColor = [0.8 0.2 0.2];
%     ax1.Box = 'off';
%     ax2.Box = 'off';
%     
%     title(['mean asymmetry in collisions for kappa=' num2str(params.K)])
% %     legend('mean asym in space','mean asym in time','Location','northwest')
%     if SAVEFIG
%         savefig(['collision asym means'])
%         saveas(h,['collision asym means' '.png'],'png')
%     end
%     
%     for i = [1 floor(Sample/2) Sample]
%         clear h
%         h = figure;
%         plot(x(1:end-1),asymmetry(i,:))
%         xlabel('x')
%         ylabel('asymmetry in collisions')
%         title(['asymmetry in collisions time-step ' num2str(i)])
%     %     legend('Raw Collisions','Smoothed Collisions')
%         if SAVEFIG
%             savefig(['collision asym time-step ' num2str(i)])
%             saveas(h,['collision asym time-step ' num2str(i) '.png'],'png')
%         end
%     end
    


end
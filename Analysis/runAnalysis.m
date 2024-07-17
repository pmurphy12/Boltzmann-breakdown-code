function [] = runAnalysis(dataArray,params,MovieTitle)
% Runs MCkernelPlot_2d_to_1d.m, beforeAfterComp.m, and ProfChange.m for
% each simulation dataset stored in dataArray. Data are assumed to be from
% repeats with parameters given in structure params.
    
%     StatsSum = cell(1,length(dataArray)+1);
    for i = 1:length(dataArray)
        [~,data1D,~,~] = MCkernelPlot_2d_to_1d(dataArray{i},params,'save_name',[MovieTitle ' ' num2str(i)],'simType',' Wedge');
        beforeAfterComp(data1D,params,[MovieTitle ' ' num2str(i)]);
        Stats = ProfChange(dataArray{i},params);
        Stats.Type = MovieTitle;
%         StatsSum{i} = Stats;
        save([MovieTitle ' ' num2str(i) ' Stats' ],'Stats')
    end
    
%     StatsSum{end} = MovieTitle;
    
    
    
    
end


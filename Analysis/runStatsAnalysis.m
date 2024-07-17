function [StatsSum] = runStatsAnalysis(kappaValueString)
% Runs MCkernelPlot_2d_to_1d.m, beforeAfterComp.m, and ProfChange.m for
% each simulation dataset stored in dataArray. Data are assumed to be from
% repeats with parameters given in structure params.
    
    
    
%     Files = dir(fullfile('./',[MovieTitle '*Wedge' '*2021.mat']));
%     Files = dir(fullfile('./',[MovieTitle '*data.mat']));
    Files = dir(fullfile('./',['*' kappaValueString '*data.mat']));
%     StatsSum = zeros(1,length(Files));
    
    for i = 1:length(Files)
        listOfVariables = who('-file',Files(i).name);
        if length(listOfVariables) ~= 2
            error(['incorrect number of variables in ' Files(i).name])
        end
        S = load(Files(i).name,listOfVariables{:});
        testString = strcmp('data',listOfVariables{1}(1:4));
        if testString
            Stats = ProfChange(S.(listOfVariables{1}),S.(listOfVariables{2}));
            f = fieldnames(S.(listOfVariables{2}));
            for ii = 1:length(f)
               Stats.(f{ii}) = S.(listOfVariables{2}).(f{ii});
            end
%             Stats.Type = MovieTitle;
            a = convertCharsToStrings(Files(i).name);
            b = convertCharsToStrings(kappaValueString);
            [~,I] = regexp(a,b);
            Stats.Type = Files(i).name(1:I);
        else
           Stats = ProfChange(S.(listOfVariables{2}),S.(listOfVariables{1}));
            f = fieldnames(S.(listOfVariables{1}));
            for ii = 1:length(f)
               Stats.(f{ii}) = S.(listOfVariables{1}).(f{ii});
            end
%             Stats.Type = MovieTitle;
            a = convertCharsToStrings(Files(i).name);
            b = convertCharsToStrings(kappaValueString);
            [~,I] = regexp(a,b);
            Stats.Type = Files(i).name(1:I);
        end
        
        StatsSum(i) = Stats;
    end
    
%     StatsSum{end} = MovieTitle;
%     save([MovieTitle ' Stats' ],'StatsSum')
    save([kappaValueString ' Stats' ],'StatsSum')
    
%     F = fieldnames(StatsSum);
% 
%     for i = 1:length(F)-1
%         A.(F{i}) = mean(cat(2,StatsSum(:).(F{i})),2);
%     end
%     save(kappaValueString ' average stats','A')
    
    
    
end


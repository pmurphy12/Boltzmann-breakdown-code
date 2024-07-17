function [dataCat,params,Stats] = runConcat_and_Stats(identificationString)
% Concatenates each simulation dataset stored in dataArray. Data are
% assumed to be from repeats with parameters given in structure params.
% Runs ProfChange.m for concatenation.

    SAVE = true;
    STATS = false;
    
    
%     Files = dir(fullfile('./',[MovieTitle '*Wedge' '*2021.mat']));
%     Files = dir(fullfile('./',[MovieTitle '*data.mat']));
    Files = dir(fullfile('./',['*' identificationString '*.mat']));
%     StatsSum = zeros(1,length(Files));
    
    % extract information for first file to gent params and length of data
    listOfVariables = who('-file',Files(1).name);
    if length(listOfVariables) ~= 3 && length(listOfVariables) ~= 2
        error(['incorrect number of variables in ' Files(1).name])
    end
    S = load(Files(1).name,listOfVariables{:});
    for ii = 1:length(listOfVariables)
        shortList1{ii} = listOfVariables{ii}(1:4);
        shortList2{ii} = listOfVariables{ii}(1:min(6,length(listOfVariables{ii})));
    end

    testString1 = strcmp('data',shortList1);
    testString2 = strcmp('params',shortList2);

    params = S.(listOfVariables{testString2});
    
    data = cell(length(Files),length(S.(listOfVariables{testString1})));
    data(1,:) = S.(listOfVariables{testString1});

    for i = 2:length(Files)
        listOfVariables = who('-file',Files(i).name);
        if length(listOfVariables) ~= 3 && length(listOfVariables) ~= 2
            error(['incorrect number of variables in ' Files(i).name])
        end
        S = load(Files(i).name,listOfVariables{:});
        for ii = 1:length(listOfVariables)
            shortList1{ii} = listOfVariables{ii}(1:4);
            shortList2{ii} = listOfVariables{ii}(1:min(6,length(listOfVariables{ii})));
        end
        
        testString1 = strcmp('data',shortList1);
        testString2 = strcmp('params',shortList2);
        
        data(i,:) = S.(listOfVariables{testString1});
        

    end

    dataCat = strucCat(data);
    
    if STATS
        Stats = ProfChange(dataCat,params);
    end

    f = fieldnames(params);
    for ii = 1:length(f)
       Stats.(f{ii}) = params.(f{ii});
    end

%         a = convertCharsToStrings(Files(i).name);
%         b = convertCharsToStrings(kappaValueString);
%         [~,I] = regexp(a,b);
%         Stats.Type = Files(i).name(1:I);

%         a = pwd;
%         [~,I] = regexp(a,'L');
%         Stats.Type = a(I:end);
%         
%         StatsSum(i) = Stats;

    if SAVE && STATS
        save(['Stats and concatenated rod info' ],'Stats','dataCat','params','-v7.3')
    elseif SAVE
        save(['Concatenated rod info' ],'dataCat','params','-v7.3')
    end
    
    
end


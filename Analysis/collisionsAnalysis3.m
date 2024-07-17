function [Stats] = collisionsAnalysis3(identificationString)
% Calculated number of collisions per time step for population and for
% individual cells. Additionally calculates mean free path over set
% interval and its standard deviation.

% Takes in data from 1D projection of agent-based simulation data.

    Title = 'full sim';
    SAVEDATA = true;

    Files = dir(fullfile('./',['*' identificationString '*.mat']));
    
    % extract information for first file to get params and length of data
    listOfVariables = who('-file',Files(1).name);
    if length(listOfVariables) ~= 3
        error(['incorrect number of variables in ' Files(1).name])
    end
    S = load(Files(1).name,listOfVariables{:});
    for ii = 1:length(listOfVariables)
        shortList1{ii} = listOfVariables{ii}(1:4);
        shortList2{ii} = listOfVariables{ii}(1:min(6,length(listOfVariables{ii})));
        shortList3{ii} = listOfVariables{ii}(1:min(10,length(listOfVariables{ii})));
    end

%     testString1 = strcmp('data',shortList1);
    testString2 = strcmp('params',shortList2);
    testString3 = strcmp('collisions',shortList3);

    params = S.(listOfVariables{testString2});
    c = S.(listOfVariables{testString3});
    %     data = S.(listOfVariables{testString1});
    
    collisions = cell(length(Files),1);
    collisions{1} = c';
    
    
    % Proceed with the rest of the files now that params structure has been
    % loaded
    for ii = 2:length(Files)
        listOfVariables = who('-file',Files(ii).name);
        if length(listOfVariables) ~= 3
            error(['incorrect number of variables in ' Files(ii).name])
        end
        S = load(Files(ii).name,listOfVariables{:});
        for iii = 1:length(listOfVariables)
            shortList1{iii} = listOfVariables{iii}(1:4);
%             shortList2{iii} = listOfVariables{iii}(1:min(6,length(listOfVariables{iii})));
            shortList3{iii} = listOfVariables{iii}(1:min(10,length(listOfVariables{iii})));
        end
        
%         testString1 = strcmp('data',shortList1);
        testString3 = strcmp('collisions',shortList3);

        c = S.(listOfVariables{testString3});
        %     data = S.(listOfVariables{testString1});

        collisions{ii} = c';
    
    end
    
    collisionsM = [collisions{:}];
    
    tSpan = 1:params.Nt;
    timing = cell(size(collisionsM,2),1);
    diffTime = timing;
    for i = 1:length(timing)
        timing{i} = tSpan(logical(collisionsM((tSpan),i)));
        diffTime{i} = diff(timing{i})*params.dt;
    end
    
    Stats.meanTime = mean([diffTime{:}]);
    Stats.stdTime = std([diffTime{:}]);
    Stats.kurtosisTime = kurtosis([diffTime{:}]);
    
    Stats.cumColls = cumsum(sum(collisionsM,2));
    Stats.cumCollsPcell = Stats.cumColls/length(Files)/params.N;
    Stats.collNum = sum(collisionsM(2:end,:),2)/length(Files);
    Stats.collNumSmoothed = movingmean(Stats.collNum,floor(5/params.dt));
    Stats.meanCollNumber = mean(sum(collisionsM,1));
    Stats.stdCollNumber = std(sum(collisionsM,1));
    Stats.kurtosisCollNumber = kurtosis(sum(collisionsM,1));

    
    if SAVEDATA
        save(['Statistics for collisions in ' Title],'Stats','collisionsM')
    end


end
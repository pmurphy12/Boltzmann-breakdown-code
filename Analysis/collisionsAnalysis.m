function [Stats] = collisionsAnalysis(identificationString,varargin)

% Takes in data from 1D projection of agent-based simulation data.

%     p = inputParser;
%     addOptional(p,'LengthX', 1, @isnumeric);
%     addOptional(p,'LengthY',1,@isnumeric);
%     addOptional(p,'save_name','',@ischar);
%     addOptional(p,'simType',' Trianglewave Test',@ischar);
%     addRequired(p,'MCdata');
%     addRequired(p,'MCparams');
    
    
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
    N = 1:params.Nt;
    
    
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
    path2length_ratio = [timing timing];
    for i = 1:length(timing)
        timing{i} = N(logical(collisionsM((tSpan),i)));
        diffTime{i} = diff(timing{i})*params.dt;
        path2length_ratio{i,1} = diffTime{i}*params.vbar/params.len;
        if length(timing{i}) > 1
            path2length_ratio{i,2} = timing{i}(2:end);
        else
            path2length_ratio{i,2} = [];
        end
    end
    
    path2length_ratioCat = [[path2length_ratio{:,1}]' [path2length_ratio{:,2}]'];
    
    for i = ceil(3*params.Nt/4):floor(4*params.Nt/4)
        p2len(i) = mean(path2length_ratioCat(path2length_ratioCat(:,2) == i,1));
    end
    
    Stats.meanTime = mean([diffTime{:}]);
    Stats.stdTime = std([diffTime{:}]);
    Stats.meanDist = Stats.meanTime*params.vbar;
    Stats.meandist_over_length = Stats.meanDist/params.len;
    
%     a = [diffTime{:}];
%     histogram(a)
%     median(a)
%     ans*params.vbar
%     ans/params.len
%     b = sort(a);
    
%     t = linspace(0,params.tend,params.Nt);

    save(['Statistics for collisions in last quarter of sim'],'Stats','timing','diffTime','collisionsM','path2length_ratioCat','p2len')


end
function [T_abs] = profileAnalysisRun(identificationString)
% Calculate table of mean-field profile properties. Meant to be run from
% the "Mean Field" folder.
% identificationString = 'kappa'

    Title = 'full sim';
    dataFileName = '1D data combined.mat';
    SAVEFIG = true;
    SAVEDATA = true;
    
%     addpath('./mean-field data')

    Folders = dir(fullfile('./',['*' identificationString '*']));
    
    varNames = {'time','peakLocation','peakAmplitude',...
        'skewness','kappa','Num','wave','yshift'};
    
    % extract information from files
%     listOfVariables = who('-file',['./mean-field data/' Files(1).name]);


    time = zeros(length(Folders),101);
    kappa = zeros(length(Folders),101);
    Num = zeros(length(Folders),101);
    peakLocation = zeros(length(Folders),101);
    peakAmplitude = zeros(length(Folders),101);
    skewness = zeros(length(Folders),101);
    wave_direction = cell(length(Folders),101);
    yshift = zeros(length(Folders),101);
    
    tableheight = 0;
    
    
    for ii = 1:length(Folders)
        load([Folders(ii).name filesep dataFileName],'data1D');
        load([Folders(ii).name filesep 'data for run 1.mat'],'params');
        data1D = cat(1,data1D{2,:});
        tt = linspace(0,params.tend,size(data1D,1));
        x = linspace(0,params.lengthGx,size(data1D,2));
        if strcmp(Folders(ii).name(48:53),'yshift')
            yshift_val = str2double(Folders(ii).name(55:56));
        else
            yshift_val = 0;
        end
            
        parfor jj = 1:length(tt)
            time(ii,jj) = tt(jj);
            kappa(ii,jj) = params.K;
            Num(ii,jj) = params.N;
            % get statistics for profiles
            [peakLocation(ii,jj),peakAmplitude(ii,jj),skewness(ii,jj)] = profileAnalysis(x,data1D(jj,:));
            wave_direction{ii,jj} = 'right-moving';
            yshift(ii,jj) = yshift_val;
                           
        end
        
        tableheight = tableheight+length(tt);
    end
    
    
    % Table generation
    
    
    T_abs = table(time(kappa ~= 0),...
                peakLocation(kappa ~= 0),...
                peakAmplitude(kappa ~= 0),...
                skewness(kappa ~= 0),...
                kappa(kappa ~= 0),...
                Num(kappa ~= 0),...
                wave_direction(kappa ~= 0),...
                yshift(kappa ~= 0),...
                'VariableNames',varNames);
            
%     T_abs = T_abs(T_abs.kappa ~= 0,:);
    
    if SAVEDATA
        save(['Agent profile features ' Title],'T_abs')
    end
    
    


end
% designed to run beforeAfterCompBins2 for selected folders


identificationString = 'scaling test';
Folders = dir(fullfile('./',['*' identificationString '*']));
dataFileName = 'Concatenated rod info.mat';

for i = 1:length(Folders)
    
    cd(Folders(i).name)
    
    load('data for run 1.mat','params');
    load(dataFileName,'dataCat');
    beforeAfterCompBins2(dataCat,params);
    
    clear params dataCat
    cd('../')
end







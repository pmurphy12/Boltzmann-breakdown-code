function [] = runMCscript_quasirefract(results_folder)

curr_fold = pwd;
addpath('./Analysis')
addpath('./Parallelogram')

if exist([results_folder],'dir') == 0
    mkdir([results_folder]);
else
    addpath(results_folder);
end

Kappa = 0.1;
N = 45700;%4570*80;%2.856e+03;%45700; %2000*2;
L = 400;
len = sqrt(2*L^2/N*Kappa);%(1/(2*sqrt(2))/4); %2*sqrt(2); 2*sqrt(5)
tend = 50; %/(sqrt(3)/2);
m = 0.1; % normal is 0.1 (l/v), quarter step is 0.4, sixteenth step is 1.6, etc. m = 1 -> dt = l/(10v)

alpha = len/sqrt(2*L^2/2856*0.1);
% alpha = 1;

Nruns = 100;

for i = 1:Nruns
    
    [params,data,collisions] = alignMC_quasirefractory(i,len,N,L,tend,m,'RaisedWave','Polar',...
        'SaveFolder',[results_folder '/'],'InitialSD',[25 25]/L,'InitialMean',[150 250]/L,...
        'yShift',false,'yShift_all',true,'yShift_supp',sqrt(alpha)*10,'NoiseLevel',0,'NoCollisions',false);
    save([results_folder '/data for run ' num2str(i)],'params','data','collisions')
    clear params data collisions
    %'yShift_supp',sqrt(alpha)*10
end

cd(results_folder)

[dataCat,params] = runConcat_and_Stats('data');
beforeAfterCompBins2(dataCat,params);
% movieCompBins2(dataCat,params);
plotWave('.');

cd(curr_fold)


end
function [] = runMyxoGen(results_folder)
    

    if exist([results_folder],'dir') == 0
        mkdir([results_folder]);
    else
        addpath(results_folder);
    end

    N = 2.856e+03;%2.856e+03;%45700;
    kappa = 0.1;%0:0.02:0.9; %0.7

    tend = 50;

    for i = 1:length(kappa)

        myxo_gen(N,kappa(i),tend,results_folder)

    end


end
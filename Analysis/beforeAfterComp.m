function [] = beforeAfterComp(data1D,params,varargin)

% Takes in data from 1D projection of agent-based simulation data.

%     p = inputParser;
%     addOptional(p,'LengthX', 1, @isnumeric);
%     addOptional(p,'LengthY',1,@isnumeric);
%     addOptional(p,'save_name','',@ischar);
%     addOptional(p,'simType',' Trianglewave Test',@ischar);
%     addRequired(p,'MCdata');
%     addRequired(p,'MCparams');


    for wave = 1:2
        L = length(data1D{wave,1});
        x = linspace(0,params.lengthGx/params.lengthGy,L);
        mean0 = circStats(x,data1D{wave,1},params.lengthGx/params.lengthGy);
        meanT = circStats(x,data1D{wave,end},params.lengthGx/params.lengthGy);
        K1 = floor((params.lengthGx/params.lengthGy/2-mean0)*L/params.lengthGx*params.lengthGy);
        K2 = floor((params.lengthGx/params.lengthGy/2-meanT)*L/params.lengthGx*params.lengthGy);
        dataShift1 = circshift(data1D{wave,1},K1);
        dataShift2 = circshift(data1D{wave,end},K2);
        
%         dataShift = circshift(data1D{wave,end},-floor(1/sqrt(2)*0.3*sqrt(2)*256));
        
        figure;
        plot(x,dataShift1*params.lengthGy,'k--')
        hold on
        plot(x,dataShift2*params.lengthGy,'b')
        xlabel('x')
        ylabel('Probability Density')
        title('Density before and after one collision')
        
        if isempty(varargin)
            savefig(['Density before and after one collision for density ' num2str(wave)])
        else
            savefig([varargin{1} ' density before and after one collision for density ' num2str(wave)])
        end
    end







end
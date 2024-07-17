function [data1D] = beforeAfterCompBins(data,params,varargin)

% Takes in data from 1D projection of agent-based simulation data.

%     p = inputParser;
%     addOptional(p,'LengthX', 1, @isnumeric);
%     addOptional(p,'LengthY',1,@isnumeric);
%     addOptional(p,'save_name','',@ischar);
%     addOptional(p,'simType',' Trianglewave Test',@ischar);
%     addRequired(p,'MCdata');
%     addRequired(p,'MCparams');

    data1D = cell(2,length(data));
    x = linspace(0,params.lengthGx,2^8);
    Nsims = length(data{1}.o)/params.N;

    angleBins = [-pi -pi/2 0 pi/2 pi];
    nBins = length(angleBins)-1;
    nCells = cell(nBins,2);
    k = 1; % index for nCells

    for i = 1:length(data)
        Pos = data{i}.pos;
        Ang = data{i}.o;
        
        for ii = 1:nBins
            pos = Pos(Ang >= angleBins(ii) & Ang < angleBins(ii+1),:);
            if i == 1 || i == length(data)
                nCells{ii,k} = size(pos,1);
                k = k+1;
            end
            data1D{ii,i} = ksdensity(pos,x,'BoundaryCorrection','reflection','Support',[0,params.lengthGx]);
        end
        
    end


    h = figure;
    for wave = 1:nBins



        L = length(data1D{wave,1});
        x = linspace(0,params.lengthGx,L);
%         mean0 = circStats(x,data1D{wave,1},params.lengthGx/params.lengthGy);
%         meanT = circStats(x,data1D{wave,end},params.lengthGx/params.lengthGy);
%         K1 = floor((params.lengthGx/params.lengthGy/2-mean0)*L/params.lengthGx*params.lengthGy);
%         K2 = floor((params.lengthGx/params.lengthGy/2-meanT)*L/params.lengthGx*params.lengthGy);
%         dataShift1 = circshift(data1D{wave,1},K1);
%         dataShift2 = circshift(data1D{wave,end},K2);
        data1 = data1D{wave,1}*nCells{wave,1}/Nsims/params.lengthGy;
        data2 = data1D{wave,end}*nCells{wave,2}/Nsims/params.lengthGy;
        
        
        plot(x,data1,'k--')
        hold on
        
        plot(x,data2,'Color',color_chooser(wave))
        

    end

    xlabel('x')
    ylabel('Number Density')
%         title('Density before and after one collision')
    title('Density before and after')
%     legend('right-moving wave initial','right-moving wave','left-moving wave initial','left-moving wave')

    if isempty(varargin)
        savefig(['Density before and after'])
        saveas(h,['Density before and after.png'],'png')
    else
        savefig([varargin{1} ' density before and after one collision for density ' num2str(wave)])
    end
%     if isempty(varargin)
%         savefig(['Density before and after one collision for density ' num2str(wave)])
%         saveas(h,['Density before and after one collision for density ' num2str(wave) '.png'],'png')
%     else
%         savefig([varargin{1} ' density before and after one collision for density ' num2str(wave)])
%     end




end
function [data1D] = beforeAfterCompBins2(data,params,varargin)

% Takes in data from 1D projection of agent-based simulation data.

%     p = inputParser;
%     addOptional(p,'LengthX', 1, @isnumeric);
%     addOptional(p,'LengthY',1,@isnumeric);
%     addOptional(p,'save_name','',@ischar);
%     addOptional(p,'simType',' Trianglewave Test',@ischar);
%     addRequired(p,'MCdata');
%     addRequired(p,'MCparams');
    Wrap = false;

    data1D = cell(2,length(data));
    x = linspace(0,params.lengthGx,2^8);
    Nsims = length(data{1}.o)/params.N;

    angleBins = [-pi -pi/2 0 pi/2 pi];
    nBins = length(angleBins)-1;
    nCells = cell(nBins);
    I = {[1 4],[2 3]}; % indicies for left- and right-moving waves respectively


    pos = cell(nBins,1);

    for i = 1:length(data)
        Pos = data{i}.pos;
        if Wrap
            Pos(:,1) = mod(Pos(:,1) + params.lengthGx/2,params.lengthGx);
        end
        Ang = data{i}.o;
        
        for ii = 1:nBins
            pos{ii} = Pos(Ang >= angleBins(ii) & Ang < angleBins(ii+1),:);
            nCells{ii} = size(pos{ii},1);
        end
        
        data1D{1,i} = ksdensity([pos{1}(:,1); pos{4}(:,1)],x,'BoundaryCorrection','reflection','Support',[0,params.lengthGx])...
            *sum([nCells{I{1}}])/Nsims/params.lengthGy;
        data1D{2,i} = ksdensity([pos{2}(:,1); pos{3}(:,1)],x,'BoundaryCorrection','reflection','Support',[0,params.lengthGx])...
            *sum([nCells{I{2}}])/Nsims/params.lengthGy;

    end
    
    save('1D data combined','data1D')


    h = figure;
    for wave = 1:2 % right-moving then left-moving waves



        L = length(data1D{wave,1});
        x = linspace(0,params.lengthGx,L);
%         mean0 = circStats(x,data1D{wave,1},params.lengthGx/params.lengthGy);
%         meanT = circStats(x,data1D{wave,end},params.lengthGx/params.lengthGy);
%         K1 = floor((params.lengthGx/params.lengthGy/2-mean0)*L/params.lengthGx*params.lengthGy);
%         K2 = floor((params.lengthGx/params.lengthGy/2-meanT)*L/params.lengthGx*params.lengthGy);
%         dataShift1 = circshift(data1D{wave,1},K1);
%         dataShift2 = circshift(data1D{wave,end},K2);
        data1 = data1D{wave,1};
        data2 = data1D{wave,end};
        
        if wave == 1
            plot(x,data1,'k--')
        else
            plot(x,data1,'k.-')
        end
        hold on
        if wave == 1
            plot(x,data2,'r')
        else
            plot(x,data2,'b')
        end

    end
    
    xlabel('x')
    ylabel('Number Density')
%         title('Density before and after one collision')
    title('Density before and after')
    legend('left-moving wave initial','left-moving wave','right-moving wave initial','right-moving wave',...
        'Location','southeast')

    if isempty(varargin)
        savefig(['Density before and after'])
        saveas(h,['Density before and after.png'],'png')
    else
        savefig([varargin{1} ' density before and after one collision for density ' num2str(wave)])
    end

    


end
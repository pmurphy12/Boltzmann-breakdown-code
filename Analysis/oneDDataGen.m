function [data1D] = oneDDataGen(data,params,varargin)

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
    nCells = cell(nBins);
    I = {[1 4],[2 3]}; % indicies for left- and right-moving waves respectively


    pos = cell(nBins,1);

    for i = 1:length(data)
        Pos = data{i}.pos;
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




end
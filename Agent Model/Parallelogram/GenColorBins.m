function [PosBinned,OrientBinned] = GenColorBins(Pos,Orient,NColorBins,varargin)
% GenColorBins creates two NColorArray length cell arrays where each cell is
% a N x 2 array of positions and orientations binned based on NColorBins 
% bins uniformly distributed in [-pi,pi]. Pos is a N x 2 array of agent 
% positions in 2d. Orient is a N x 2 array where the 
% i-th row is an orientation vector for agent i = 1...N.
%
% TO DO: Add varargin to also bin based on input vector BinLocations.
    
    dims = size(Orient);
    N = dims(1); %Number of agents.

    if ~isempty(varargin)
        Bins = varargin{1};
        NColorBins = length(Bins)-1;
    elseif isempty(varargin)
        Bins = linspace(-pi,pi,NColorBins+1);
    else
        error('Exceeded the maximum of 4 inputs.')
    end
    
    if length(dims) ~= 2
        error('Incorrect number of dimensions for data. Needs to be N x M.')
    end
    CBins = false(N,NColorBins);

    PosBinned = cell(NColorBins,1);
    OrientBinned = cell(NColorBins,1);
    Ang = atan2(Orient(:,2),Orient(:,1));

    for i = 1:NColorBins

        CBins(:,i) = (Ang >= Bins(i)).*(Ang < Bins(i+1));
        OrientBinned{i} = Orient(CBins(:,i),:);
        PosBinned{i} = Pos(CBins(:,i),:);

    end


    




end


function [m] = findMode(x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    [h,I] = histcounts(x);
    [~,J] = max(h);
    m = I(J);
    
end


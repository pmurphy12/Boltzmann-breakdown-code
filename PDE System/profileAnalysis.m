function [peakLocation,peakAmplitude,skewness,standDev] = profileAnalysis(x,y)
%   Calculate interquartile range, skewness, and change in front and back
%   of support above EPS for profiles.
    
%     L2norm = zeros(length(ANG),length(data));
%     RENORMALIZE y?
    [peakAmplitude,Imax] = max(y);
    peakLocation = x(Imax);
    
    m1 = moment(y,1);
    m2 = moment(y,2);
    m3 = moment(y,3);
    skewness = (m3-3*m1*m2+2*m1^3)/(m2-m1^2)^1.5;
    standDev = sqrt(m2 - m1^2);
    
    
    
    
end


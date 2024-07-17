function [cMean,cVar,cSTD,cDisp] = circStatsSample(posData,L)
%circVar calculated the circular variacne, standard deviation, and
%dispertion for a dataset of samples (posData) on a 1D domain  of length L
%with periodic boundary conditions.
    
    %Converte to points on the circle
    z = exp(1i*posData*2*pi/L);
    m1 = mean(z);
    m2 = mean(z.^2);
    
    % Calculate mean on domain [0,L].
    cMean = L/(2*pi)*mod(angle(m1),2*pi);
    
    % Calculate Variance
    R = abs(m1);
    cVar = 1-R;
    
    %Calculate Standard Deviation
    cSTD = sqrt(-2*log(R));
    
    %Calculate circular dispersion
    R2 = abs(m2);
    cDisp = (1-R2)/(2*R^2);
    
    
    
end


function [cMean,cVar,cSTD,cDisp] = circStats(xi,distData,L)
%circVar calculated the circular variacne, standard deviation, and
%dispertion for a dataset of distrbution values (distData) on a 1D domain
%of length L with periodic boundary conditions. xi is a vector of the
%spatial locations corresponding to distData.
    
    %Converte to points on the circle
    z = exp(1i*xi(1:end-1)*2*pi/L);
    m1 = 2*pi*mean(distData(1:end-1).*z);
    m2 = 2*pi*mean(distData(1:end-1).*(z.^2));
%     m1 = 2*pi/(length(xi))*sum(distData(1:end).*z(1:end));
%     m2 = 2*pi/(length(xi))*sum(distData(1:end).*(z(1:end).^2));
    
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


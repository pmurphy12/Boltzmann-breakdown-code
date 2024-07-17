function [cMean,cVar,cSTD] = eucStats(xi,distData)
%circVar calculated the circular variacne, standard deviation, and
%dispertion for a dataset of distrbution values (distData) on a 1D domain
%of length L with periodic boundary conditions. xi is a vector of the
%spatial locations corresponding to distData.
    
    distData = distData/((xi(2)-xi(1))*sum(distData(1:end-1)));
    
    m1 = (xi(2)-xi(1))*sum(distData(1:end-1).*xi(1:end-1));
%     m2 = (xi(2)-xi(1))*sum(distData(1:end-1).*xi(1:end-1).^2);
    m2 = (xi(2)-xi(1))*sum(distData(1:end-1).*(xi(1:end-1)-m1).^2);
%     m1 = 2*pi/(length(xi))*sum(distData(1:end).*z(1:end));
%     m2 = 2*pi/(length(xi))*sum(distData(1:end).*(z(1:end).^2));
    
    % Calculate mean on domain [0,L].
    cMean = m1;
    
    % Calculate Variance
    cVar = m2;
    
    %Calculate Standard Deviation
    cSTD = sqrt(cVar);
    
    
    
end


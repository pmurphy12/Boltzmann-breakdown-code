
function [yprime, tmes] = rhsm(t,y,dx,c1,c2,N,K,tmes,tmesout)

yprime = zeros(2*N,1);


r1 = y(1:N);
r2 = y(N+1:2*N);


if abs(t - tmes - tmesout) < 0.0001
    fprintf('t= %f\n',t);
    tmes = tmes + tmesout;
end

f1   =  c1 * fluxf(r1)   .* fluxg(r2,c2);
f2   = -c1 * fluxf(r2)   .* fluxg(r1,c2);

% Low-order scheme

kp1 = [2:N 1]';
km1 = [N 1:N-1]';

% this is F_j+1/2
F1 = (f1(kp1) + f1)/2;
F2 = (f2(kp1) + f2)/2;

du1 = r1(kp1) - r1;
du2 = r2(kp1) - r2;

% 
F1 = F1 - K/2 * du1;
F2 = F2 - K/2 * du2;

%if ilim == 1
  % Limiter
%end

% Final result
yp1 = -(F1 - F1(km1))/dx;
yp2 = -(F2 - F2(km1))/dx;


yprime(1:N) = yp1;
yprime(N+1:2*N) = yp2;





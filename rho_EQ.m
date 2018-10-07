function [Rv,sum_c] = rho_EQ(Vs,D,V)
%% here D = D^2, sigma = sigma^2;
Rv = V;
vT = 1; vR = 0;  indp = find(V>vR);
sqrtD = sqrt(D);
intovT = dawson((vT-Vs)/sqrtD)*exp((vT-Vs)^2/D);
% intovT,
% pause;
intovSD = dawson(-Vs/sqrtD)*exp(Vs^2/D);
% Vs^2/D, Vs, D,
%% compute R with V>vR case:
Rv(indp) = - dawson((V(indp)-Vs)/sqrtD) ...
    + exp(-(V(indp) - Vs).^2/D)* intovT;

if (indp(1)>1)
    Rv(1:indp(1)-1) = exp(-(V(1:indp(1)-1) - Vs).^2/D)* (-intovSD + intovT);
end
indp = find(V<-2/3);
Rv(1:indp) = 0;
sum_c = (V(2)-V(1))*sum(Rv);
Rv = Rv/sum_c;
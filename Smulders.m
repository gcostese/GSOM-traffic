function v=Smulders(r,I)
%Computation of the speed according to a given spacing (r)
%and a given driver attribute (I) = class of the vehicle
%with the generalized Smulders model
% I==1 for cars
% I==2 for trucks

rho=1/r; %passage à la densité (inverse du spacing)

Rjam=   1/5.*(I==1) + 1/15.*(I==2);   %densité maximale (en veh/m)
Vmax=    34.*(I==1) + 25.*(I==2);     %vitesse libre (en m/s)
Rcrit= 1/30.*(I==1) + 1/35.*(I==2);   %densité critique (en veh/m)
Vcrit=   22.*(I==1) + 20.*(I==2);     %vitesse critique (en m/s)

beta=(Vmax-Vcrit)/Rcrit;
w=(Rcrit * Vcrit)/(Rjam-Rcrit);


v= max(0 + Vmax.*(rho==0) , ...
    (Vmax-beta*rho).*(rho <= Rcrit) ...
    + w*(Rjam/rho-1).*(rho > Rcrit) );
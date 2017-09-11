function mu=partialSmulders(r,I)
%Calcul de la valeur de la dérivée (mu) du diagramme fondamental à la
%Smulders, pour la valeur de spacing (r) et pour l'attribut
%conducteur (I)

rho=1/r;

Rjam=   1/5.*(I==1) + 1/15.*(I==2);   %densité maximale (en veh/m)
Vmax=    34.*(I==1) + 25.*(I==2);     %vitesse libre (en m/s)
Rcrit= 1/30.*(I==1) + 1/35.*(I==2);   %densité critique (en veh/m)
Vcrit=   22.*(I==1) + 20.*(I==2);     %vitesse critique (en m/s)

beta=(Vmax-Vcrit)/Rcrit;
w=(Rcrit * Vcrit)/(Rjam-Rcrit);

mu= (beta*rho^2).*(rho <= Rcrit) ...
    + w*Rjam.*(rho > Rcrit) ;
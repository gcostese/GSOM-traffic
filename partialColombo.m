function mu=partialColombo(r,I)
%Calcul de la valeur de la dérivée (mu) du diagramme fondamental à la
%Colombo 1-phase, pour la valeur de spacing (r) et pour l'attribut
%conducteur (I)

rho=1/r; %passage à la densité (inverse du spacing)
R=1/5; %densité maximale (en veh/m)
V=25; %vitesse libre (en m/s)
Rhat=1/30; %densité critique (en veh/m)
Vhat=22; %vitesse critique (en m/s)
beta=(V-Vhat)/Rhat; %paramètre de pente de la fonction vitesse en fluide
qstar=1; %paramètre

A=V+qstar/R-I;
B=beta-I/R;
rhocrit=1/(2*B)*(A-sqrt(A^2-4*qstar*B));%densité critique

if rho <= rhocrit
    mu=beta*rho^2;
else
    mu=qstar+I*rho^2/R;
end
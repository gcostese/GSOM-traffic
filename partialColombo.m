function mu=partialColombo(r,I)
%Calcul de la valeur de la d�riv�e (mu) du diagramme fondamental � la
%Colombo 1-phase, pour la valeur de spacing (r) et pour l'attribut
%conducteur (I)

rho=1/r; %passage � la densit� (inverse du spacing)
R=1/5; %densit� maximale (en veh/m)
V=25; %vitesse libre (en m/s)
Rhat=1/30; %densit� critique (en veh/m)
Vhat=22; %vitesse critique (en m/s)
beta=(V-Vhat)/Rhat; %param�tre de pente de la fonction vitesse en fluide
qstar=1; %param�tre

A=V+qstar/R-I;
B=beta-I/R;
rhocrit=1/(2*B)*(A-sqrt(A^2-4*qstar*B));%densit� critique

if rho <= rhocrit
    mu=beta*rho^2;
else
    mu=qstar+I*rho^2/R;
end
function v=Colombo(r,I)
%Calcul de la valeur de la vitesse pour un spacing (r) donné et un attribut
%véhicule (I) fixé, selon le modèle de Colombo 1-phase

rho=1/r; %passage à la densité (inverse du spacing)
R=1/5; %densité maximale (en veh/m)
V=25; %vitesse libre (en m/s)
Rhat=1/30; %densité critique (en veh/m)
Vhat=22; %vitesse critique (en m/s)
beta=(V-Vhat)/Rhat; %paramètre de pente de la fonction vitesse en fluide
qstar=1; %paramètre
%Imax=qstar/R; %indice maximal

A=V+qstar/R-I;
B=beta-I/R;
rhocrit=1/(2*B)*(A-sqrt(A^2-4*qstar*B));%densité critique

if rho <= rhocrit %partie fluide
    v=V-beta*rho;
else              %partie congestionnée
    v=(I+qstar/rho)*(1-rho/R);
end
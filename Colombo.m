function v=Colombo(r,I)
%Calcul de la valeur de la vitesse pour un spacing (r) donn� et un attribut
%v�hicule (I) fix�, selon le mod�le de Colombo 1-phase

rho=1/r; %passage � la densit� (inverse du spacing)
R=1/5; %densit� maximale (en veh/m)
V=25; %vitesse libre (en m/s)
Rhat=1/30; %densit� critique (en veh/m)
Vhat=22; %vitesse critique (en m/s)
beta=(V-Vhat)/Rhat; %param�tre de pente de la fonction vitesse en fluide
qstar=1; %param�tre
%Imax=qstar/R; %indice maximal

A=V+qstar/R-I;
B=beta-I/R;
rhocrit=1/(2*B)*(A-sqrt(A^2-4*qstar*B));%densit� critique

if rho <= rhocrit %partie fluide
    v=V-beta*rho;
else              %partie congestionn�e
    v=(I+qstar/rho)*(1-rho/R);
end
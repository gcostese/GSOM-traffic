function r=invSmulders(v,I)
%fuonction permettant de calculer la valeur inverse (r) par le diagramme
%fondamental à la Smulders, en connaissant la valeur de vitesse (v)
%et l'attribut conducteur (I)

Rjam=   1/5.*(I==1) + 1/15.*(I==2);   %densité maximale (en veh/m)
Vmax=    34.*(I==1) + 25.*(I==2);     %vitesse libre (en m/s)
Rcrit= 1/30.*(I==1) + 1/35.*(I==2);   %densité critique (en veh/m)
Vcrit=   22.*(I==1) + 20.*(I==2);     %vitesse critique (en m/s)

% if v>Vmax
%     disp('Error: speed exceeds maximal speed!')
% end

beta=(Vmax-Vcrit)/Rcrit;
w=(Rcrit * Vcrit)/(Rjam-Rcrit);


rho=(Vmax-v)/beta.*(v >= Vcrit)...
    + (w*Rjam)/(w+v).*(v < Vcrit);

r=1/rho;
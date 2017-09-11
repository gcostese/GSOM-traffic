function L=Lagrangian_Smulders(mu,I)
%fonction qui permet de donner le "coût" (ou Lagrangien)
%d'un point de la trajectoire, localement de pente v
%le Lagrangien est la transformée de Legendre de l'Hamiltonien qui est
%défini par la fonction Smulders(r,I)
%L(v,I)=sup_r {Smulders(r,I)-rv}

L=0;

Rjam=   1/5.*(I==1) + 1/15.*(I==2);   %densité maximale (en veh/m)
Vmax=    34.*(I==1) + 25.*(I==2);     %vitesse libre (en m/s)
Rcrit= 1/30.*(I==1) + 1/35.*(I==2);   %densité critique (en veh/m)
Vcrit=   22.*(I==1) + 20.*(I==2);     %vitesse critique (en m/s)

beta=(Vmax-Vcrit)/Rcrit;
w=(Rcrit * Vcrit)/(Rjam-Rcrit);

rcrit=1/Rcrit;

v1=sqrt(beta/mu);

g1=@(r) Vmax-beta/r-mu*r;
g2=@(r) (w*Rjam-mu)*r-w;

L= Smulders(inf,I).*(mu==0) ...
    + max( ...
    g1(v1).*(rcrit<=v1) ...
    + g1(rcrit).*(rcrit>v1) , ...
    g2(rcrit) ).*(mu>0) ;
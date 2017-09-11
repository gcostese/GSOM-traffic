function r=invColombo(v,I)
%fuonction permettant de calculer la valeur inverse (r) par le diagramme
%ondamental à la Colombo 1-phase, en connaissant la valeur de vitesse (v)
%et l'attribut conducteur (I)

R=1/5; %densité maximale (en veh/m)
V=25; %vitesse libre (en m/s)
Rhat=1/30; %densité critique (en veh/m)
Vhat=22; %vitesse critique (en m/s)
beta=(V-Vhat)/Rhat; %paramètre de pente de la fonction vitesse en fluide
qstar=1; %paramètre

A=V+qstar/R-I;
B=beta-I/R;
rhocrit=1/(2*B)*(A-sqrt(A^2-4*qstar*B));%densité critique
vcrit=Colombo(1/rhocrit,I);

if v >= vcrit %fluide
    r=beta/(V-v);
else %congestionné
    a=qstar*R;
    b=(I-v)*R-qstar;
    c=-I;
    delta=b^2-4*a*c;
    r=(-b+sqrt(delta))/(2*a);
end
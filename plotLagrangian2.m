function plotLagrangian2(vmax)
%fonction permettant de tracer sur un même diagramme, le Lagrangien de W
%(DF speed spacing) pour les attributs I=1, I=2 et I=3
%entre les valeurs 0 et vmax
g=figure;
hold on

pas=0.05;
for I=0:5
    %initialisation de la matrice
    M=zeros(1,fix(vmax/pas));
    p=1;
    for v=0:pas:vmax
        %calcul de la valeur de la vitesse W pour chacun des attributs
        M(1,p)=Lagrangian2(v,I);
        p=p+1;
    end
    plot((0:pas:vmax),M(1,:),'b','LineWidth',1)
end
xlabel('u (1/s)','Fontsize',12)
ylabel('Cost (m/s^2)','Fontsize',12)
title('Lagrangian M(N,u,t)','Fontsize',14)
hold off

saveas(g,'lagrangian.eps','epsc')
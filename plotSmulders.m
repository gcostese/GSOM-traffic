function plotSmulders(rmin,rmax)
%fonction permettant de tracer sur un même diagramme, la fonction de
%vitesse pour les attributs I, entre les valeurs rmin et rmax
z=figure;
hold on

pas=0.05;
for I=[1 2]
    %initialisation de la matrice
    M=zeros(1,fix((rmax-rmin)/pas));
    p=1;
    for r=rmin:pas:rmax
        %calcul de la valeur de la vitesse W pour chacun des attributs
        M(1,p)=Smulders(r,I);
        p=p+1;
    end
    plot((rmin:pas:rmax),M(1,:),'b','LineWidth',1)
end
xlabel('Spacing r (m)','Fontsize',12)
ylabel('Speed W (m/s)','Fontsize',12)
title('Fundamental diagram W(r,I)','Fontsize',14)
hold off

saveas(z,'fundamentaldiag-Smulders.eps','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tracé du Diagramme Fondamental débit-densité associé
q=figure;
hold on

rmax=20*rmax;
pas=0.05;
for I=[1 2]
    %initialisation de la matrice
    M=zeros(1,fix((rmax-rmin)/pas)+1);
    p=1;
    for r=rmin:pas:rmax
        %calcul de la valeur de la vitesse W pour chacun des attributs
        M(1,p)=3600*Smulders(r,I)/r;
        p=p+1;
    end
    plot(1000./(rmin:pas:rmax),M(1,:),'r','LineWidth',1)
end
xlabel('Density \rho (veh/km)','Fontsize',12)
ylabel('Flow F (veh/h)','Fontsize',12)
title('Fundamental diagram F(\rho,I)','Fontsize',14)
hold off

saveas(q,'fundamentaldiag2-Smulders.eps','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tracé du Diagramme Fondamental vitesse-densité associé
q=figure;
hold on

rmax=20*rmax;
pas=0.05;
for I=[1 2]
    %initialisation de la matrice
    M=zeros(1,fix((rmax-rmin)/pas)+1);
    p=1;
    for r=rmin:pas:rmax
        %calcul de la valeur de la vitesse W pour chacun des attributs
        M(1,p)=Smulders(r,I);
        p=p+1;
    end
    plot(1000./(rmin:pas:rmax),M(1,:),'m','LineWidth',1)
end
xlabel('Density \rho (veh/km)','Fontsize',12)
ylabel('Speed W (m/s)','Fontsize',12)
title('Fundamental diagram V(\rho,I)','Fontsize',14)
hold off

saveas(q,'fundamentaldiag3-Smulders.eps','epsc')
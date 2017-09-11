function [h,g]=plotSolution(X,DeltaN2,Deltat2,N0,t0,N1,N2,t1,t2,A,B,Neuler,Teuler,test)
%fonction permettant de tracer la solution finale entre les indices N1 et
%N2 et entre les temps t1 et t2

%To dock AUTOMATICALLY the figures
%set(0,'DefaultFigureWindowStyle','docked') 
P12=fix((N2-N1)/DeltaN2)+1;
P22=fix((t2-t1)/Deltat2)+1;

N=linspace(N1,N2,P12);
t=linspace(t1,t2,P22);
[N,t]=meshgrid(N,t);

%Deltat2=Deltat/10
%DeltaN2=DeltaN/10
i1=fix((t1-t0)/Deltat2)+1;
i2=fix((t2-t0)/Deltat2)+1;
j1=fix((N1-N0)/DeltaN2)+1;
j2=fix((N2-N0)/DeltaN2)+1;
Xp=X(i1:i2,j1:j2);
%Xp=X(i1:i2,j2:-1:j1);

h=figure;
mesh(N,t,Xp,'FaceColor','interp','FaceLighting','phong');
camlight right
colormap winter
%colorbar
xlabel('Vehicle label','FontSize',16)
ylabel('Time (s)','FontSize',16)
zlabel('Position (m)','FontSize',16)

hold on
if ~isempty(B)==1 && strcmp(test,'yes2')==1
    for x0=B
        %Tracer la courbe des véhicules cumulés sur le plot 3D
        i3=find(Neuler==N1);
        i4=find(Neuler==N2);
        Neuler=Neuler(i3:i4);
        Teuler=Teuler(i3:i4);
        s=length(Neuler);
        Xeuler=(x0+5).*ones(1,s);
        plot3(Neuler,Teuler,Xeuler,'-r','LineWidth',2);
    end
end
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tracer la valeur des vitesses sur le pavé (N,t)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure supplémentaire= trajectoires des véhicules (labels are integers)
g=figure;
hold on
V=(N1:N2);
if strcmp(test,'yes')==1
    for n=V(~ismember(V,A)) %the labels of A are excluded
        plot((t1:Deltat2:t2)',X(i1:i2,fix((n-N0)/DeltaN2)+1));
    end
    
    %Plot in red the modified trajectories (in case of Lagrangian data
    %assimilation)
    for n=A
        if n>=N1 && n<=N2
            plot((t1:Deltat2:t2)',X(i1:i2,fix((n-N0)/DeltaN2)+1),'r');
        end
    end
else
    for n=V
        plot((t1:Deltat2:t2)',X(i1:i2,fix((n-N0)/DeltaN2)+1));
    end
end
ylabel('Location (m)','FontSize',16)
xlabel('Time (s)','FontSize',16)
title('Vehicles trajectories','FontSize',18)
hold off
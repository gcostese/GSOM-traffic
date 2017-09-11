function [Xn0,vn0,Rn0]=CinternLagr2(t0,tmax,Deltat,x0,i0)
%fonction qui permet de d�terminer la position du v�hicule (Xn0) 
%d'indice n_0, pour tout temps t in [t0,tmax]
%vn0 est la vitesse instantann�e du v�hicule d'indice n_0
%(attention � �tre coh�rent avec la vitesse maximale donn�e pour la
%fonction W et admissible pour le v�hicule n_0

P2=fix((tmax-t0)/Deltat)+1;

Xn0=zeros(P2,1);
%--------------------------------------------------------------------------
% vmax=25;
% vn0=vmax.*rand(P2,1);
% --> A remplacer par: trajectoire d'un v�hicule en d�c�l�ration
v1=18;
v2=10;
vn0=Speed(v1,v2,t0,tmax,Deltat)';
%--------------------------------------------------------------------------
Rn0=zeros(P2,1);
Xn0(1)=x0;
Rn0(1)=invColombo(vn0(1),i0);
for p=2:P2
    Xn0(p)=Xn0(p-1)+vn0(p-1)*Deltat;
    Rn0(p)=invColombo(vn0(p),i0);
end
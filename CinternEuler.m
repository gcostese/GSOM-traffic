function [N,T,V,R,F]=CinternEuler(T0,Pini,Pfin,DeltaN,N0,tmax,I)
%fonction qui permet de d�terminer l'�quation de N(t) telle que nous avons
%       X(N(t),t)=x0   (condition "interne" eul�rienne)
%La courbe N(t) est consid�r�e comme affine par morceaux i.e. la pente
%\dot{N}(t) qui est homog�ne � un flux, est constante par morceaux
%Cette pente doit respecter la condition de compatibilit� � savoir que la
%pente est inf�rieure au flux maximal donn� pour l'attribut v�hicule I(N,t)

% ------- OUTPUT ----------
%T = matrix of the times such that for each t_p \in T, we have
% N(t_p)=N_p for p \in [Pini,Pfin]
%V = matrix of the speeds on the trajectory N(t)
%R = matrix of the spacings transported on the trajectory N(t)

%Parameters
rmin=5;
rmax=50;

%Instantiation
T=zeros((Pfin-Pini+1),1);
V=zeros((Pfin-Pini),1);
F=zeros((Pfin-Pini),1);

%Spacing
R=rmin+(rmax-rmin).*rand((Pfin-Pini+1),1); %en m

T(1,1)=T0;
for p=Pini:Pfin
    %Label
    N=(p-1)*DeltaN+N0;
    
    %Attribute
    In=I(N);
    
    %Speed
    r=R(p-Pini+1,1);
    v=Colombo(r,In); %en m/s
    V(p-Pini+1,1)=v;
    
    %Flow
    f=1/r*v; %en 1/s
    F(p-Pini+1,1)=f;
    
    %Time
    T(p-Pini+2,1)=T(p-Pini+1,1)+DeltaN/f;
    
    Nfin=N;
    if T(p-Pini+2)>tmax
        break
    end
end

Nini=N0+(Pini-1)*DeltaN;
N=(Nini:Nfin);

%Resizing
s=length(N);

T=T(1:s,1);
R=R(1:s,1);
V=V(1:s,1);
F=F(1:s,1);

% %------- V�rification graphique -------
% %1) Condition interne
% figure;
% Nini=N0+(Pini-1)*DeltaN;
% Nfin=N0+Pfin*DeltaN;
% plot(Nini:Nfin,T)
% 
% %2) Diagramme Fondamental D�bit-Densit�
% figure;
% plot(1./R,F,'.')
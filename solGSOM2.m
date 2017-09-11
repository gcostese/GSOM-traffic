function [X,Xmod,Xmod2]=solGSOM2(N0,Nmax,t0,tmax)
%sch�ma de r�solution propos� par JP Lebacque afin de r�duire le co�t de
%calcul par rapport � la m�thode plus �l�gante, exacte mais "it�rative"
%Approche matricielle (plus performante sous Matlab)

%Calcul parall�le
matlabpool open
tic

%% Param�tres
%[t0 tmax]: plage de temps consid�r�
%[N0 Nmax]: plage de v�hicules consid�r�s
DeltaN=1; %pas en indice (ici v�hicules consid�r�s 1 par 1)
P1=fix((Nmax-N0)/DeltaN); %nombre d'intervalles de discr�tisation en indice
Deltat=1; %pas de temps (si conditions limites et mixtes)
P2=fix((tmax-t0)/Deltat); %nombre d'intervalles de discr�tisation en temps
P=P1+P2; %nombre d'intervalles de conditions (ini + bound) � consid�rer
%P1+P2+1 points (P1+1 en indices et P2+1 en temps, -1 commun en (t0,N0))

%Colombo(N,r,t)=diagramme fondamental vitesse-interdistance (r)

%% Conditions initiales (pour t=t_0)
%Xt0 positions des v�hicules � l'instant initial X(t_0,N) 
%pour tout N \in [N_0,N_max]
%Anciennement [Rt0,It0,Xt0]=CI(N0,Nmax);
rmin=5; %interdistance minimale en m�tres
rmax=50; %interdistance "maximale" en m�tres
Xt0=zeros(P1+1,1);
Rt0=(rmax-rmin).*rand(P1+1,1)+rmin;
It0=5.*rand(P1+1,1);
Xt0(1)=0;
for p=2:P1+1
    Xt0(p)=Xt0(p-1)-Rt0(p-1)*DeltaN;
end

%I(N,t)=attribut li� aux conducteurs qui r�sout �quation de transport
%|\dot{I}(N,t)=\phi(I)
%|I(t_0)=I_0
%Anciennement I=@(n) dynI(It0(n));
I=@(n) It0(fix((n-N0)/DeltaN)+1);

%% Conditions aux bords (pour N=N_0)
%Xn0 positions du v�hicule N_0 pour tout temps t \in [t_0,t_max]
%Anciennement [Xn0,vn0]=CB(t0,tmax);
%          et Rn0=@(t) invColombo(vn0(t),I(N0));
Xn0=zeros(P2+1,1);
%--------------------------------------------------------------------------
% vmax=25;
% vn0=vmax.*rand(P2+1,1);
% --> A remplacer par: trajectoire d'un v�hicule en acc�l�ration
v1=15;
v2=20;
vn0=Speed(v1,v2,t0,tmax,Deltat)';
%--------------------------------------------------------------------------
Rn0=zeros(P2+1,1);
Xn0(1)=0;
Rn0(1)=invColombo(vn0(1),It0(1));
for p=2:P2+1
    Xn0(p)=Xn0(p-1)+vn0(p-1)*Deltat;
    Rn0(p)=invColombo(vn0(p),It0(1));
end
V0=vn0(1);


yon='y';
%yon=input('Visualisation (''y'' or ''n'') : ');
if strcmp(yon,'y')==1
    %trac� du DF entre spacing 5 et 200 (en m�tres)
    plotColombo(5,100)
    %trac� du Lagrangien entre 0 et mu=3
    plotLagrangian2(2.5);
    %repr�sentation graphique des conditions initiales
    plotCI2(N0,Nmax,DeltaN,Rt0,I,Xt0)
    %repr�sentation graphique des conditions au bord
    plotCB2(t0,tmax,Deltat,Xn0,vn0,Rn0)
    close all
end

%% Conditions internes lagrangiennes (pour un ensemble d'indices N=N_i)
%Xni positions du v�hicule N_i pour tout temps t \in [t_0,t_max]
%M�me cas que ci-dessous mais simplement translat�!
A=MatAlea(N0,Nmax,5);
%A=input('Indices des v�hicules (conditions internes): ');

%% Conditions internes Euleriennes (pour X=x_0)
%Cumulated Vehicle Count at a fixed (Eulerian) point for any time t in
%the adequate set of time include in [t_0,t_max]
% x0=min(Xt0);
% xmax=max(Xn0);
% B=MatAlea(x0,xmax,0.5);
B=0;
% B=input('Position de la boucle de comptage (conditions eul�riennes): ');


%% Calcul des solutions partielles
%Ins�rer une grille plus fine (pour interpolation)
DeltaN2=DeltaN/10;
Deltat2=Deltat/10;
%pour tous les multiples (N0+k*10*DeltaN2), on est sur un N=N_q
%idem en temps

%Initialisation matrice finale de la solution globale
P12=fix((Nmax-N0)/DeltaN2)+1;
P22=fix((tmax-t0)/Deltat2)+1;
%temps verticalement, indices horizontalement
X=zeros(P22,P12);
X(:,:)=inf;

%Positions initiales = 1ere ligne de la matrice finale
col=1;
for n=N0:DeltaN2:Nmax
    p=fix((n-N0)/DeltaN)+1;
    X(1,col)=Xt0(p)-(n-(p-1)*DeltaN)*Rt0(p);
    col=col+1;
end
%Positions du v�hicule N0 (condition au bord)
%                   = 1ere colonne de la matrice finale
ligne=1;
for t=t0:Deltat2:tmax
    q=fix((t-t0)/Deltat)+1;
    X(ligne,1)=Xn0(q)+(t-(q-1)*Deltat)*vn0(q);
    ligne=ligne+1;
end

%Slice the matrix X to avoid unnecessary communication overhead
Mini=X(1,:);
Mbound=X(:,1);

%parfor p=1:P
for p=1:P
    
    %Initialisation de la matrice des solutions partielles �mises par la
    %condition "p"
    %temps verticalement, indices horizontalement
    Xp=zeros(P22,P12);
    Xp(:,:)=inf;
    if p<=P1
        %% Cas des conditions initiales
        %______________________________________________________________________
        %Initialisation de la solution partielle pour t0
        Np0=N0+(p-1)*DeltaN; %indice de d�but de l'intervalle
        rp0=Rt0(p); %valeur de l'interdistance sur l'intervalle
        Xp(1,:)=Mini;

        %Calcul pour les conditions initiales (p=1:P1)
        %calcul de 2 ou 3 (si rar�faction) caract�ristiques aux points {N_q} 
        %pour tout q>=p et tel que N_q<=Nmax
        %indexer les caract�ristiques pour $i$ \in \{1,2\} ou \{1,2,3\}
        %pour chaque intervalle [p,p+1], d�terminer les temps T{i}(p,q) pour 
        %lesquels la caract�ristique (i) atteint un indice N_q (q>=p) o� I est 
        %sens� �voluer
        %______________________________________________________________________
        %D�terminer le nombre de caract�ristiques (=savoir s'il y a
        %rar�faction ou non pour le premier point de chaque intervalle)
        %1. d�terminer les caract�ristiques de l'�tat de l'intervalle pr�c�dent
        %(� "gauche" d'o� la d�nomination)
        if p==1
            vgauche=V0;
        else
            rgauche=Rt0(p-1);
            Igauche=I(Np0-DeltaN2);
            vgauche=Colombo(rgauche,Igauche);
        end
        %D�terminer le spacing projet� sur le DF de la bande consid�r�e
        if vgauche<Colombo(inf,I(Np0))
            rstar=invColombo(vgauche,I(Np0));
        else
            rstar=inf;
        end
        if rstar>rp0
            nbcaract=3; %cas de l'onde de rar�faction
        else
            nbcaract=2; %cas sans onde de rar�faction
        end
        %disp([num2str(nbcaract),' caracteristiques'])
        %__________________________________________________________________
        %Initialisation des matrices auxiliaires
        %matrice des temps de passage de chaque caract�ristique
        %aux fronti�res des bandes [n_p,n_{p+1}]
        T=zeros(nbcaract,1);
        %matrice des pentes \partial_r V(r,I)
        mu=zeros(nbcaract,1);
        %matrice des spacings r(n,t)
        R=zeros(nbcaract,1);
        %matrice des co�ts (L pour Lagrangien) point � point
        L=zeros(nbcaract,1);
        %matrice des valeurs positions en chaque point de passage
        Y=zeros(nbcaract,1);
        %__________________________________________________________________
        %Initialisation des temps de passage et des spacings pour la 
        %premi�re cellule [n_p,n_{p+1}]
        %Correspond au cas $q=p$
        T(1:nbcaract,1)=t0;
        %__________________________________________________________________
        %Calcul des caract�ristiques, de la valeur de X sur chacune d'elle,
        %puis calcul dans le domaine d'influence (calcul exact entre les
        %caract�ristiques 1 et 2 MAIS interpolation pour l'�ventail de
        %rar�faction i.e. entre les caract�ristiques 2 et 3)
        Nold=Np0;
        Told=t0;
        q=1;
        while Nold<Nmax && Told<tmax     
            if q==1
                Nnew=Np0;
                %initialisation des spacings transport�s par la 
                %caract�ristique $i$ dans l'intervalle [n_p,n_{p+1}]
                j=round((Np0-N0)/DeltaN2)+1;
                Y(1,1)=Mini(1,j+floor(DeltaN/DeltaN2));
                Y(2:nbcaract,1)=Mini(1,j);
                for i=1:nbcaract
                    if i==3
                        R(3,1)=rstar; %le spacing transport� par la caract.
                                      %3 est diff�rent
                        if rstar==inf
                            T(3,2)=tmax; %attention, cas d�g�n�r�
                            disp('Error: infinite spacing value!')
                            %break
                        end
                    else
                        R(i,1)=rp0;
                    end
                    mu(i,1)=partialColombo(R(i,1),I(Np0));
                    T(i,2)=t0+DeltaN/mu(i,1);
                    if i==1
                        T(i,2)=t0;
                    end
                    L(i,1)=Lagrangian2(mu(i,1),I(Np0));
                    Y(i,q+1)=(T(i,2)-t0)*L(i,1)+Y(i,q);
                end
            else %q>=2
                Nnew=Nold+DeltaN;
                %Calcul par morceaux de l'�quation de la caract�ristique
                %ainsi que de la solution Xp le long des caract�ristiques
                for i=1:nbcaract
                    Rold=R(i,q-1);
                    if Colombo(Rold,I(Nold))<Colombo(inf,I(Nnew))
                        %on se projete sur le nouveau DF si possible
                        Rnew=invColombo(Colombo(Rold,I(Nold)),I(Nnew));
                    else
                        %sinon on est dans le cas critique et le spacing est infini
                        %Rnew=inf;
                        disp('Error: infinite spacing value!')
                        %break
                    end
                    %Actualisation des variables
                    R(i,q)=Rnew; %valeur du spacing
                    mu(i,q)=partialColombo(Rnew,I(Nnew)); %pente
                    T(i,q+1)=T(i,q)+DeltaN/mu(i,q); %temps travers�e

                    %Formule de Lax-Hopf
                    %d�terminer la position X{p}(Nq,T{i}(p,q))
                    %=position du point Nq au temps d�termin� ci-dessus T{i}(p,q)
                    L(i,q)=Lagrangian2(mu(i,q),I(Nnew));
                    %calcul de la position au temps de passage en n_{q+1}
                    Y(i,q+1)=(T(i,q+1)-T(i,q))*L(i,q)+Y(i,q);
                end
            end
            Told=T(1,q+1);
            Nold=Nnew;
            
            q=q+1; %actualisation de l'indice de la "bande"
            fprintf('.')
        end %boucle "while" (domaine de computation)
        qmax=q-1;
        
        %__________________________________________________________________
        %Etape 2: Calculer les positions aux points de la grille secondaire
        for q=1:qmax
            Nq=Np0+(q-1)*DeltaN; %indice de base
            %a=[Nq T(2,q)]';
            vo=Y(2,q);
            %disp(vo)
            spacing=R(1,q);
            %disp(spacing)
            speed=Colombo(spacing,I(Nq));

            for j=1:10
                n=Nq+j*DeltaN2; %sous-indice
                ind=(p+q-2)*DeltaN/DeltaN2+j+1;
                
                %Calcul uniquement dans la bande de caract�ristiques
                pente=1/mu(1,q);
                Tmin=min(max(0,T(1,q+1)-pente*(10-j)*DeltaN2),tmax);
                % r�gle le cas du premier intervalle (q=1)
                Tmax=min(T(2,q)+pente*j*DeltaN2,tmax);
                t1=T(2,q)+(n-Nq)/mu(2,q);
                %disp([n Tmin Tmax t1])
                tpmin=ceil((Tmin-t0)/Deltat2)+1;
                tpmax=floor((Tmax-t0)/Deltat2)+1;
                for tps=tpmin:1:tpmax
                    %Calcul exact pour les points situ�s dans la bande
                    %des caract�ristiques
                    t=t0+(tps-1)*Deltat2;
                    vx=vo+L(2,q)*(t1-T(2,q))+speed*(t-t1);
                    %disp(vx)
                    Xp(tps,ind)=vx;
                end

                if nbcaract==3
                    %disp('rarefaction')
                    %Interpolation de la valeur � chaque point de la sous-grille 
                    %� partir des valeurs aux points de grille
                    %M�thode: interpolation barycentrique � partir des 3 plus 
                    %proches voisins, en ne consid�rerant uniquement que les points
                    %dans un rectangle pertinent
                    
                    pentemin=1/mu(2,q);
                    pentemax=1/mu(3,q);
                    Tmin=min(T(2,q)+pentemin*j*DeltaN2,tmax);
                    Tmax=min(T(3,q)+pentemax*j*DeltaN2,tmax);
                    tpmin=ceil((Tmin-t0)/Deltat2)+1;
                    tpmax=floor((Tmax-t0)/Deltat2)+1;
                    if q==1 %cas particulier du premier intervalle
                        for tps=tpmin:1:tpmax
                            t=t0+(tps-1)*Deltat2;
                            x=[n t]';
                            vx=inf;

                            %1 seul triangle (avec rar.) � consid�rer
                            a=[Nq T(2,q)]';
                            b=[Nq+DeltaN T(2,q+1)]';
                            c=[Nq+DeltaN T(3,q+1)]';
                            va=Y(2,q);
                            vb=Y(2,q+1);
                            vc=Y(3,q+1);
                            vx=min(vx,interpol(x,a,b,c,va,vb,vc));
                            
                            %OldX=Xp(tps,ind);
                            Xp(tps,ind)=vx; %min(OldX,vx);
                        end
                    else %tous les cas diff�rents de q==1
                        for tps=tpmin:1:tpmax
                            t=t0+(tps-1)*Deltat2;
                            x=[n t]';
                            vx=inf;
                            
                            for cas=1:2 %2 triangles si 3 ondes
                                a=[Nq T(3,q)]';
                                b=[Nq+DeltaN T(2,q+1)]';
                                va=Y(3,q);
                                vb=Y(2,q+1);
                                if mod(cas,2)==1 %triangles impairs
                                    c=[Nq T(2,q)]';
                                    vc=Y(2,q);
                                else %mod(cas,2)==0, triangles pairs
                                    c=[Nq+DeltaN T(3,q+1)]';
                                    vc=Y(3,q+1);
                                end
                                vx=min(vx,interpol(x,a,b,c,va,vb,vc));
                            end %boucle "for" sur les cas de triangles
                            
                            %OldX=Xp(tps,ind);
                            Xp(tps,ind)=vx; %min(OldX,vx);
                        end %boucle"for" sur les temps
                        
                    end %boucle "if" sur la valeur de q
                    
                end %boucle "if" pour tester si rar�faction
                
            end %boucle "for" sur l'indice j pour Nq+j*DeltaN2
        
        end %boucle "for" sur l'indice q de la bande consid�r�e

        %%Test logique pour savoir si les ondes se recouvrent bien
        %if p>=2 && T{p}(nbcaract,2)<T{p-1}(1,2)
        %    disp(['Pas de recouvrement pour p=',num2str(p)]);
        %    break
        %end

        %%%%%%%%%%%%%%%%%%%%%%%%% TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Xp=Xp(1:P22,1:P12); %matrix dimension must agree
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp([' End for p= ',num2str(p)]);
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    elseif p>P1 && p<=P
        %% Cas des conditions aux bords
        
        %Calcul sym�trique pour les conditions au bord � N=N0 (t=1:P2)
        pp=p-P1;
        %__________________________________________________________________
        %Initialisation de la solution partielle pour N0
        tp0=t0+(pp-1)*Deltat; %indice de d�but de l'intervalle
        rp0=Rn0(pp); %valeur de l'interdistance sur l'intervalle
        Xp(:,1)=Mbound;

        Np0=N0; %indice de base auquel est accroch� la donn�e lagrangienne

        %Calcul pour les conditions initiales (pp=1:P2)
        %calcul de 2 ou 3 (si rar�faction) caract�ristiques aux points {N_q} 
        %pour tout q>=p et tel que N_q<=Nmax
        %indexer les caract�ristiques pour $i$ \in \{1,2\} ou \{1,2,3\}
        %pour chaque intervalle [p,p+1], d�terminer les temps T{i}(p,q) pour 
        %lesquels la caract�ristique (i) atteint un indice N_q (q>=p) o� I est 
        %sens� �voluer
        %__________________________________________________________________
        %D�terminer le nombre de caract�ristiques (=savoir s'il y a
        %rar�faction ou non pour le premier point de chaque intervalle)
        %1. d�terminer les caract�ristiques de l'�tat de l'intervalle pr�c�dent
        %(� "gauche" d'o� la d�nomination)
        if pp==1
            rstar=Rt0(1);
        else %pp>2
            rstar=Rn0(pp-1);
        end
        %D�terminer le nombre de caract�ristiques � calculer en fonction de r*
        if rstar<rp0
            nbcaract=3; %cas de l'onde de rar�faction
        else
            nbcaract=2; %cas sans onde de rar�faction
        end
        %disp([num2str(nbcaract),' caracteristiques'])% pour p=',num2str(p)]);
        %______________________________________________________________________
        %matrice des valeurs positions en chaque point de passage
        Y=zeros(nbcaract,1);
        %matrice des temps de passage de chaque 
        %caract�ristique aux fronti�res des bandes [n_p,n_p+1]
        T=zeros(nbcaract,1);
        %matrice des pentes \partial_r V(r,I)
        mu=zeros(nbcaract,1);
        %matrice des spacings r(n,t)
        R=zeros(nbcaract,1);
        %matrice des co�ts (L pour Lagrangien) point � point
        L=zeros(nbcaract,1);
        %______________________________________________________________________
        %Initialisation des temps de passage et des spacings pour la premi�re
        %cellule [N_0,N_1]
        %Correspond au cas $q=1$
        T(1,1)=tp0+Deltat;
        T(2:nbcaract,1)=tp0;
        %______________________________________________________________________
        %Calcul des caract�ristiques, de la valeur de X sur chacune d'elle,
        %puis calcul dans le domaine d'influence (calcul exact entre les
        %caract�ristiques 1 et 2 MAIS interpolation pour l'�ventail de
        %rar�faction i.e. entre les caract�ristiques 2 et 3)
        Nold=Np0;
        Told=tp0;
        q=1;
        while Nold<Nmax && Told<tmax
            if q==1
                Nnew=Np0;
                %initialisation des spacings transport�s par la caract�ristique
                % $i$ dans l'intervalle [N_0,N_1]
                i=fix((tp0-t0)/Deltat2)+1;
                Y(1,q)=Mbound(i+floor(Deltat/Deltat2),1);
                Y(2:nbcaract,q)=Mbound(i,1);
                for i=1:nbcaract
                    if i==3
                        R(3,1)=rstar; %le spacing transport� par la caract. 3
                    else
                        R(i,1)=rp0;
                    end
                    mu(i,1)=partialColombo(R(i,1),I(Np0));
                    T(i,2)=T(i,1)+DeltaN/mu(i,1);
                    L(i,1)=Lagrangian2(mu(i,1),I(Np0));
                    Y(i,q+1)=(T(i,2)-T(i,1))*L(i,1)+Y(i,q);
                end
            else %q>=2
                Nnew=Np0+(q-1)*DeltaN;
                %Calcul par morceaux de l'�quation de la caract�ristique
                %ainsi que de la solution Xp le long des caract�ristiques
                for i=1:nbcaract
                    Rold=R(i,q-1);
                    if Colombo(Rold,I(Nold))<Colombo(inf,I(Nnew))
                        %on se projete sur le nouveau DF si possible
                        Rnew=invColombo(Colombo(Rold,I(Nold)),I(Nnew));
                    else
                        %sinon on est dans le cas critique et le spacing 
                        %est infini
                        %Rnew=inf;
                        disp('Error: infinite spacing value!')
                        %break
                    end
                    %Actualisation des variables
                    R(i,q)=Rnew; %valeur du spacing
                    mu(i,q)=partialColombo(Rnew,I(Nnew)); %pente
                    T(i,q+1)=T(i,q)+DeltaN/mu(i,q); %temps travers�e

                    %Formule de Lax-Hopf
                    %d�terminer la position X{p}(Nq,T{i}(p,q))
                    %=position du point Nq au temps d�termin� ci-dessus 
                    % T{i}(p,q)
                    L(i,q)=Lagrangian2(mu(i,q),I(Nnew));
                    %calcul de la position au temps de passage en n_{q+1}
                    Y(i,q+1)=(T(i,q+1)-T(i,q))*L(i,q)+Y(i,q);
                end
            end
            Told=T(nbcaract,q+1);
            Nold=Nnew;
            
            q=q+1; %actualisation de l'indice de la "bande"
            fprintf('.')
        end
        qmax=q-1;
        
        %__________________________________________________________________
        %Etape 2: Calculer les positions aux points de la grille secondaire
        for q=1:qmax;
            Nq=Np0+(q-1)*DeltaN; %indice de base
            %a=[Nq T(2,q)]';
                
            vo=Y(2,q);
            spacing=R(2,q);
            speed=Colombo(spacing,I(Nq));
            %disp(spacing)
            
            for j=1:10
                n=Nq+j*DeltaN2; %sous-indice
                ind=(q-1)*DeltaN/DeltaN2+j+1;
                
                %Calcul uniquement dans la bande de caract�ristiques
                pente=1/mu(1,q);
                Tmin=min(T(2,q)+pente*j*DeltaN2,tmax);
                Tmax=min(T(1,q)+pente*j*DeltaN2,tmax);
                tpmin=ceil((Tmin-t0)/Deltat2)+1;
                tpmax=floor((Tmax-t0)/Deltat2)+1;
                %disp([Tmin Tmax])
                for tps=tpmin:1:tpmax
                    %Calcul exact pour les points situ�s dans la bande
                    %des caract�ristiques
                    t=t0+(tps-1)*Deltat2;
                    t1=T(2,q)+(n-Nq)/mu(2,q);
                    vx=vo+L(2,q)*(t1-T(2,q))+speed*(t-t1);
                    %disp(vx)
                    
                    %OldX=Xp(tps,ind);
                    Xp(tps,ind)=vx; %min(OldX,vx);
                end

                if nbcaract==3
                    %Interpolation de la valeur � chaque point de la sous-grille 
                    %� partir des valeurs aux points de grille
                    %M�thode: interpolation barycentrique � partir des 3 plus 
                    %proches voisins, en ne consid�rerant uniquement que les points
                    %dans un rectangle pertinent
                    
                    pentemin=1/mu(3,q);
                    pentemax=1/mu(2,q);
                    Tmin=min(T(3,q)+pentemin*j*DeltaN2,tmax);
                    Tmax=min(T(2,q)+pentemax*j*DeltaN2,tmax);
                    tpmin=ceil((Tmin-t0)/Deltat2)+1;
                    tpmax=floor((Tmax-t0)/Deltat2)+1;
                    %disp([Tmin Tmax])
                    if q==1 %cas particulier du premier intervalle
                        a=[Nq T(2,q)]';
                        b=[Nq+DeltaN T(3,q+1)]';
                        c=[Nq+DeltaN T(2,q+1)]';
                        va=Y(2,q);
                        vb=Y(3,q+1);
                        vc=Y(2,q+1);
                        %disp([va vb vc])
                        for tps=tpmin:1:tpmax
                            t=t0+(tps-1)*Deltat2;
                            x=[n t]';
                            vx=interpol(x,a,b,c,va,vb,vc);
                            %disp(vx)
                            
                            %OldX=Xp(tps,ind);
                            Xp(tps,ind)=vx; %min(OldX,vx);
                        end
                    else %tous les cas q>1
                        for tps=tpmin:1:tpmax
                            t=t0+(tps-1)*Deltat2;
                            x=[n t]';
                            vx=inf; %permet d'�viter de r��crire la valeur du min
                            
                            for cas=1:2 %2 triangles
                                a=[Nq T(2,q)]';
                                va=Y(2,q);
                                b=[Nq+DeltaN T(3,q+1)]';
                                vb=Y(3,q+1);
                                if mod(cas,2)==1 %triangle impair
                                    c=[Nq+DeltaN T(2,q+1)]';
                                    vc=Y(2,q+1);
                                else %mod(cas,2)==0, triangle pair
                                    c=[Nq T(3,q)]';
                                    vc=Y(3,q);
                                end
                                vx=min(vx,interpol(x,a,b,c,va,vb,vc));
                            end
                            
                            %OldX=Xp(tps,ind);
                            Xp(tps,ind)=vx; %min(OldX,vx);
                        end
                    end
                end
            end
        end

        %%Test logique pour savoir si les ondes se recouvrent bien
        %if pp>=2 && T{p}(nbcaract,2)>T{p-1}(1,2)
        %    disp(['Pas de recouvrement pour pp=',num2str(pp)]);
        %    break
        %end

        %%%%%%%%%%%%%%%%%%%%%%%% TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Xp=Xp(1:P22,1:P12); %matrix dimension must agree
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp([' End for p= ',num2str(p)]);
        
    end
    
    %% Inf-morphism property
    X=min(X,Xp); %comparaison des �l�ments 2 par 2
end

%% Cas des conditions internes lagrangiennes
Xmod=X;
for Np0=A
    if Np0<=Nmax
        Xmod=LagrData2(Xmod,Np0,N0,Nmax,t0,tmax,DeltaN,Deltat,I,Rt0);
    end
end

%% Cas des conditions internes eul�riennes
Xmod2=Xmod;
Neuler=[];
Teuler=[];
for x0=B
    if x0>=min(Xt0) && x0<=max(Xn0)
        [Xmod2,Neuler,Teuler]=EulerData(Xmod2,x0,N0,Nmax,t0,tmax,DeltaN,Deltat,I);
    end
end

%% Repr�sentation graphique de la solution globale
set(0,'DefaultFigureWindowStyle','docked')

N1=N0;
N2=Nmax;
t1=t0;
t2=tmax;
test='no';
[h,g]=plotSolution(X,DeltaN2,Deltat2,N0,t0,N1,N2,t1,t2,A,B,Neuler,Teuler,test);
%Enregistrement des sorties graphiques
saveas(h,'solution.eps','epsc')
saveas(g,'trajectories.eps','epsc')
f=plotSpeed(X,Deltat2,t0,tmax,N0,Nmax);
saveas(f,'speed.eps','epsc')

if ~isempty(A)
    test='yes';
    [h,g]=plotSolution(Xmod,DeltaN2,Deltat2,N0,t0,N1,N2,t1,t2,A,B,Neuler,Teuler,test);
    %Enregistrement des sorties graphiques
    saveas(h,'solutionMod.eps','epsc')
    saveas(g,'trajectoriesMod.eps','epsc')
    f=plotSpeed(Xmod,Deltat2,t0,tmax,N0,Nmax);
    saveas(f,'speedMod.eps','epsc')
end


if ~isempty(B)
    test='no';
    [h,g]=plotSolution(Xmod2,DeltaN2,Deltat2,N0,t0,N1,N2,t1,t2,A,B,Neuler,Teuler,test);
    %Enregistrement des sorties graphiques
    saveas(h,'solutionMod2.eps','epsc')
    saveas(g,'trajectoriesMod2.eps','epsc')
    f=plotSpeed(Xmod2,Deltat2,t0,tmax,N0,Nmax);
    saveas(f,'speedMod2.eps','epsc')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[Order,InftY]=DetectAbnorm(X);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
toc
matlabpool close
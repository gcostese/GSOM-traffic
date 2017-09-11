function [Z,N_int,T_int]=EulerData(X,x0,N0,Nmax,t0,tmax,DeltaN,Deltat,I)
%EulerData: function to determine the position value on the computational
%grid for Eulerian data (cumulative vehicle count at a specified position)
%in case of random piecewise constant value for the flow along the 
%trajectory and Colombo 1-phase Fundamental Diagram

%Initialization of the matrix
Z=X;

%% PARAMETRES
%x0=position à laquelle est calculée la courbe de véhicules cumulés (CVC)
% = donnée eulerienne
DeltaN2=DeltaN/10;
Deltat2=Deltat/10;
P12=fix((Nmax-N0)/DeltaN2)+1;
P22=fix((tmax-t0)/Deltat2)+1;
P1=fix((Nmax-N0)/DeltaN)+1;

%Diagramme fondamental
R=1/5; %densité maximale (en veh/m)
V=25; %vitesse libre (en m/s)
Rhat=1/30; %densité critique (en veh/m)
Vhat=22; %vitesse critique (en m/s)
beta=(V-Vhat)/Rhat; %paramètre de pente de la fonction vitesse en fluide
qstar=1; %paramètre

A=@(I) V+qstar/R-I;
B=@(I) beta-I/R;
Rcrit=@(I) 1/(1/(2*B(I))*(A(I)-sqrt(A(I)^2-4*qstar*B(I)))); %spacing crit

%% Définition de la courbe de véhicules cumulés (CVC)
%Indices min et max entre lesquels on va enregistrer la CVC
Pini=1+floor(0.2*(P1-1)*rand(1));
Pfin=P1-floor(0.2*(P1-1)*rand(1));

Pini2=(Pini-1)*DeltaN/DeltaN2+1;

[row,~]=find(sort(X(:,Pini2))>=x0);
tm=min(row); %indice de la ligne
T0=t0+(tm-1)*Deltat2; %temps réel
%Definition of the internal boundary condition (Eulerian condition)
[N_int,T_int,V_int,R_int,F_int]=CinternEuler(T0,Pini,Pfin,DeltaN,N0,tmax,I);
%T_int = matrix of the times such that for each t_p \in T_intern
% N(t_p)=N_p for p \in [Pini,Pfin]
%V_int = matrix of the speeds on the trajectory N(t)
%R_int = matrix of the spacings transported on the trajectory N(t)
%F_int = matrix of the corresponding flows

%Resizing
s=length(N_int);
Pfin=min(Pini+s-2,Pfin);

%% Computation of the domain of influence
%parfor p=Pini:Pfin %if parallel computing
for p=Pini:Pfin
    %Initialisation de la matrice des solutions partielles émises par la
    %condition "p" (ici condition interne eulérienne)
    %temps verticalement, indices horizontalement
    Xp=zeros(P22,P12);
    Xp(:,:)=inf;
    
    %______________________________________________________________________
    %Initialisation de la solution partielle
    Np0=N0+(p-1)*DeltaN; %indice de début de l'intervalle
    rp0=R_int(p-Pini+1); %valeur de l'interdistance sur l'intervalle
    tp0=T_int(p-Pini+1); %temps de départ des premières ondes en Np0
    
    %______________________________________________________________________
    Ip0=I(Np0);
    rcrit=Rcrit(Ip0);
    
    %______________________________________________________________________
    %Déterminer le nombre de caractéristiques (=savoir s'il y a
    %raréfaction ou non pour le premier point de chaque intervalle)
    %1. déterminer les caractéristiques de l'état de l'intervalle précédent
    %(à "gauche" d'où la dénomination)
    if p==Pini
        vgauche=V_int(1);
    else
        rgauche=R_int(p-Pini);
        Igauche=I(Np0-DeltaN2);
        vgauche=Colombo(rgauche,Igauche);
    end
    %Déterminer le spacing projeté sur le DF de la bande considérée
    if vgauche<Colombo(inf,I(Np0))
        rstar=invColombo(vgauche,I(Np0));
    else
        rstar=inf;
    end
    if (rstar>rp0 && rp0>rcrit) || (rstar<rp0 && rp0<rcrit)
        nbcaract=3; %cas de l'onde de raréfaction
    else
        nbcaract=2; %cas sans onde de raréfaction
    end
    %disp([num2str(nbcaract),' caracteristiques'])
    %__________________________________________________________________
    %Initialisation des matrices auxiliaires
    %matrice des temps de passage de chaque caractéristique
    %aux frontières des bandes [n_p,n_{p+1}]
    T=zeros(nbcaract,1);
    %matrice des pentes \partial_r V(r,I)
    mu=zeros(nbcaract,1);
    %matrice des spacings r(n,t)
    R=zeros(nbcaract,1);
    %matrice des coûts (L pour Lagrangien) point à point
    L=zeros(nbcaract,1);
    %matrice des valeurs positions en chaque point de passage
    Y=zeros(nbcaract,1);
    %__________________________________________________________________
    %Initialisation des temps de passage et des spacings pour la 
    %première cellule [n_p,n_{p+1}]
    %Correspond au cas $q=p$
    T(1:nbcaract,1)=tp0;
    %__________________________________________________________________
    %Calcul des caractéristiques, de la valeur de X sur chacune d'elle,
    %puis calcul dans le domaine d'influence (calcul exact entre les
    %caractéristiques 1 et 2 MAIS interpolation pour l'éventail de
    %raréfaction i.e. entre les caractéristiques 2 et 3)
    Nold=Np0;
    Told=tp0;
    q=1;
    while Nold<Nmax && Told<tmax     
        if q==1
            Nnew=Np0;
            %initialisation des spacings transportés par la 
            %caractéristique $i$ dans l'intervalle [n_p,n_{p+1}]
            Y(1:nbcaract,1)=x0;
            for i=1:nbcaract
                if i==3
                    R(3,1)=rstar; %le spacing transporté par la caract.
                                  %3 est différent
                    if rstar==inf
                        T(3,2)=tmax; %attention, cas dégénéré
                        disp('Error: infinite spacing value!')
                        %break
                    end
                else
                    R(i,1)=rp0;
                end
                mu(i,1)=partialColombo(R(i,1),I(Np0));
                T(i,2)=tp0+DeltaN/mu(i,1);
                if i==1
                    T(i,2)=T_int(p-Pini+2);
                    Y(i,2)=x0;
                else
                    L(i,1)=Lagrangian2(mu(i,1),I(Np0));
                    Y(i,q+1)=(T(i,2)-tp0)*L(i,1)+Y(i,q);
                end
            end
        else %q>=2
            Nnew=Nold+DeltaN;
            %Calcul par morceaux de l'équation de la caractéristique
            %ainsi que de la solution Xp le long des caractéristiques
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
                T(i,q+1)=T(i,q)+DeltaN/mu(i,q); %temps traversée

                %Formule de Lax-Hopf
                %déterminer la position X{p}(Nq,T{i}(p,q))
                %=position du point Nq au temps déterminé ci-dessus T{i}(p,q)
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
    pente0=1/F_int(p-Pini+1);
    
    %Différencier selon si :
    % solution sous-critique (propagation par dessus => 
    %                           méthode type conditions initiales)
    % ou solution sur-critique (propagation par dessous =>
    %                           méthode type conditoons aux bords)
    
    if rp0>rcrit %Solution sous-critique (propagation par dessus)
        
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

                %Calcul uniquement dans la bande de caractéristiques
                pente=1/mu(1,q);
                Tmin=min(max(0,max(T(1,2)-pente0*(10-j)*DeltaN2,T(1,q+1)-pente*(10-j)*DeltaN2)),tmax);
                % règle le cas du premier intervalle (q=1)
                Tmax=min(T(2,q)+pente*j*DeltaN2,tmax);
                t1=T(2,q)+(n-Nq)/mu(2,q);
                %disp([n Tmin Tmax t1])
                tpmin=ceil((Tmin-t0)/Deltat2)+1;
                tpmax=floor((Tmax-t0)/Deltat2)+1;
                
                %disp([tpmin,tpmax])
                
                for tps=tpmin:1:tpmax
                    %Calcul exact pour les points situés dans la bande
                    %des caractéristiques
                    t=t0+(tps-1)*Deltat2;
                    vx=vo+L(2,q)*(t1-T(2,q))+speed*(t-t1);
                    %disp(vx)
                    Xp(tps,ind)=vx;
                end

                if nbcaract==3
                    %disp('rarefaction')
                    %Interpolation de la valeur à chaque point de la sous-grille 
                    %à partir des valeurs aux points de grille
                    %Méthode: interpolation barycentrique à partir des 3 plus 
                    %proches voisins, en ne considérerant uniquement que les points
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

                            %1 seul triangle (avec rar.) à considérer
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
                    else %tous les cas différents de q==1
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

                end %boucle "if" pour tester si raréfaction

            end %boucle "for" sur l'indice j pour Nq+j*DeltaN2

        end %boucle "for" sur l'indice q de la bande considérée
    
    else %Solution sur-critique (Propagation par en-dessous)
        
        for q=1:qmax;
            Nq=Np0+(q-1)*DeltaN; %indice de base
            %a=[Nq T(2,q)]';
                
            vo=Y(2,q);
            spacing=R(2,q);
            speed=Colombo(spacing,I(Nq));
            %disp(spacing)
            
            for j=1:10
                n=Nq+j*DeltaN2; %sous-indice
                ind=(p+q-2)*DeltaN/DeltaN2+j+1;
                
                %Calcul uniquement dans la bande de caractéristiques
                pente=1/mu(1,q);
                Tmin=min(max(0,T(2,q)+pente*j*DeltaN2),tmax);
                Tmax=min(max(T(1,1)+pente0*j*DeltaN2,T(1,q)+pente*j*DeltaN2),tmax);
                tpmin=ceil((Tmin-t0)/Deltat2)+1;
                tpmax=floor((Tmax-t0)/Deltat2)+1;
                
                %disp([tpmin,tpmax])
                
                for tps=tpmin:1:tpmax
                    %Calcul exact pour les points situés dans la bande
                    %des caractéristiques
                    t=t0+(tps-1)*Deltat2;
                    t1=T(2,q)+(n-Nq)/mu(2,q);
                    vx=vo+L(2,q)*(t1-T(2,q))+speed*(t-t1);
                    %disp(vx)
                    
                    %OldX=Xp(tps,ind);
                    Xp(tps,ind)=vx; %min(OldX,vx);
                end

                if nbcaract==3
                    %Interpolation de la valeur à chaque point de la sous-grille 
                    %à partir des valeurs aux points de grille
                    %Méthode: interpolation barycentrique à partir des 3 plus 
                    %proches voisins, en ne considérerant uniquement que les points
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
                            vx=inf; %permet d'éviter de réécrire la valeur du min
                            
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
                        end %boucle"for" sur les temps

                    end %boucle "if" sur la valeur de q

                end %boucle "if" pour tester si raréfaction

            end %boucle "for" sur l'indice j pour Nq+j*DeltaN2

        end %boucle "for" sur l'indice q de la bande considérée
    
    end %boucle "if" pour tester si solution sur ou sous-critique
    
    Xp=Xp(1:P22,1:P12); %matrix dimension must agree
    
    disp([' End for p = ',num2str(p)]);
    
    %Actualisation de la solution globale (inf-morphism)
    Z=min(Z,Xp); %comparaison des éléments 2 par 2
end
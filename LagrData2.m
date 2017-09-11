function Z=LagrData2(X,Np0,N0,Nmax,t0,tmax,DeltaN,Deltat,I,Rt0)
%LagrData2: function to determine the position value on the computational
%grid for Lagrangian data (trajectory of a specified vehicle)
%in case of random value for the speed along the trajectory
%and Colombo Fundamental Diagram


%PARAMETRES
%Np0=indice de base auquel est accrochée la donnée lagrangienne
DeltaN2=DeltaN/10;
Deltat2=Deltat/10;
P12=fix((Nmax-N0)/DeltaN2)+1;
P22=fix((tmax-t0)/Deltat2)+1;
P2=fix((tmax-t0)/Deltat)+1;

%Initialization of the matrix
Z=X;

%Définition de la trajectoire
K=fix((Np0-N0)/DeltaN2)+1;
%Indices min et max entre lesquels on va enregistrer la nouvelle
%trajectoire
Pini=1+floor(0.2*(P2-1)*rand(1));
Pfin=P2-floor(0.2*(P2-1)*rand(1));
%Temps minimum et maximum en valeur absolue entre lesquels on modifie la
%trajectoire (en tant que donnée lagrangienne)
tm=t0+(Pini-1)*Deltat;
tpl=t0+Pfin*Deltat;
ligne0=fix((tm-t0)/Deltat2)+1;
x0=Z(ligne0,K);
i0=I(Np0);
[Xn0,vn0,Rn0]=CinternLagr2(tm,tpl,Deltat,x0,i0);
%Xn0 positions du véhicule N_0 pour tout temps t \in [tm,tpl]

%Représentation graphique des conditions internes (Lagrangiennes)
plotCIntern(Np0,tm,tpl,Deltat,Xn0,vn0,Rn0)

%Positions du véhicule Np0 (condition interne lagrangienne)
%                   = K-ieme colonne de la matrice finale
ligne=ligne0;
for t=tm:Deltat2:tpl
    q=fix((t-tm)/Deltat)+1;
    Z(ligne,K)=Xn0(q)+(t-tm-(q-1)*Deltat)*vn0(q);
    ligne=ligne+1;
end

Z=Z(1:P22,1:P12);

%Slice the matrix X to avoid unnecessary communication overhead
Mbound=Z(:,K);

%parfor p=Pini:Pfin %if parallel computing
for p=Pini:Pfin
    %Initialisation de la matrice des solutions partielles émises par la
    %condition "p" (ici condition au bord)
    %temps verticalement, indices horizontalement
    Xp=zeros(P22,P12);
    Xp(:,:)=inf;
    %______________________________________________________________________
    %Initialisation de la solution partielle pour N0
    tp0=t0+(p-1)*Deltat; %indice de début de l'intervalle
    rp0=Rn0(p-Pini+1); %valeur de l'interdistance sur l'intervalle
    Xp(:,K)=Mbound(1:P22,1);
    
    %Calcul pour les conditions initiales (p=1:P2)
    %calcul de 2 ou 3 (si raréfaction) caractéristiques aux points {N_q} 
    %pour tout q>=p et tel que N_q<=Nmax
    %indexer les caractéristiques pour $i$ \in \{1,2\} ou \{1,2,3\}
    %pour chaque intervalle [p,p+1], déterminer les temps T{i}(p,q) pour 
    %lesquels la caractéristique (i) atteint un indice N_q (q>=p) où I est 
    %sensé évoluer
    %______________________________________________________________________
    %Déterminer le nombre de caractéristiques (=savoir s'il y a
    %raréfaction ou non pour le premier point de chaque intervalle)
    %1. déterminer les caractéristiques de l'état de l'intervalle précédent
    %(à "gauche" d'où la dénomination)
    if p==Pini
        if Pini==1
            k=fix((Np0-N0)/DeltaN)+1;
            rstar=Rt0(max(1,k-1));
        else
            rstar=Rn0(1);
        end
    else
        rstar=Rn0(p-Pini);
    end
    %Déterminer le nombre de caractéristiques à calculer en fonction de r*
    if rstar<rp0
        nbcaract=3; %cas de l'onde de raréfaction
    else
        nbcaract=2; %cas sans onde de raréfaction
    end
    %disp([num2str(nbcaract),' caracteristiques pour p=',num2str(p)]);
    %______________________________________________________________________
    %matrice des valeurs positions en chaque point de passage
    Y=zeros(nbcaract,1);
    %matrice des temps de passage de chaque 
    %caractéristique aux frontières des bandes [n_p,n_p+1]
    T=zeros(nbcaract,1);
    %matrice des pentes \partial_r V(r,I)
    mu=zeros(nbcaract,1);
    %matrice des spacings r(n,t)
    R=zeros(nbcaract,1);
    %matrice des coûts (L pour Lagrangien) point à point
    L=zeros(nbcaract,1);
    %______________________________________________________________________
    %Initialisation des temps de passage et des spacings pour la première
    %cellule [Np0,Np0+DeltaN]
    %Correspond au cas $q=p$
    T(1,1)=tp0+Deltat;
    T(2:nbcaract,1)=tp0;
    %______________________________________________________________________
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
            %initialisation des spacings transportés par la caractéristique
            % $i$ dans l'intervalle [N_0,N_1]
            i=fix((tp0-t0)/Deltat2)+1;
            Y(1,q)=Mbound(i+floor(Deltat/Deltat2),1);
            Y(2:nbcaract,q)=Mbound(i,1);
            for i=1:nbcaract
                if i==3
                    R(3,1)=rstar; %le spacing transporté par la caract. 3 
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
            %Calcul par morceaux de l'équation de la caractéristique
            %ainsi que de la solution Xp le long des caractéristiques
            for i=1:nbcaract
                Rold=R(i,q-1);
                if Colombo(Rold,I(Nold))<Colombo(inf,I(Nnew))
                    %on se projete sur le nouveau DF si possible
                    Rnew=invColombo(Colombo(Rold,I(Nold)),I(Nnew));
                else
                    %sinon on est dans le cas critique et le spacing est 
                    %infini
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
                %=position du point Nq au temps déterminé ci-dessus 
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
    for q=1:qmax
        Nq=Np0+(q-1)*DeltaN; %indice de base
        %a=[Nq T(2,q)]';
                
        vo=Y(2,q);
        spacing=R(2,q);
        speed=Colombo(spacing,I(Nq));
        %disp(spacing)

        for j=1:10
            n=Nq+j*DeltaN2; %sous-indice
            ind=K+(q-1)*DeltaN/DeltaN2+j;

            %Calcul uniquement dans la bande de caractéristiques
            pente=1/mu(1,q);
            Tmin=min(T(2,q)+pente*j*DeltaN2,tmax);
            Tmax=min(T(1,q)+pente*j*DeltaN2,tmax);
            tpmin=ceil((Tmin-t0)/Deltat2)+1;
            tpmax=floor((Tmax-t0)/Deltat2)+1;
            %disp([Tmin Tmax])
            for tps=tpmin:tpmax
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
                    end
                end
            end
        end
    end
    Xp=Xp(1:P22,1:P12); %matrix dimensions must agree
    
    disp([' End for p = ',num2str(p)]);
    
    %Actualisation de la solution globale (inf-morphism)
    Z=min(Z,Xp); %comparaison des éléments 2 par 2
end
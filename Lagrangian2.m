function L=Lagrangian2(mu,I)
%fonction qui permet de donner le "coût" (ou Lagrangien)
%d'un point de la trajectoire, localement de pente v
%le Lagrangien est la transformée de Legendre de l'Hamiltonien qui est
%défini par la fonction Colombo(r,I)
%L(v,I)=sup_r {Colombo(r,I)-rv}

L=0;

R=1/5; %densité maximale (en veh/m)
rmin=1/R; %spacing minimal (en m)
V=25; %vitesse libre (en m/s)
Rhat=1/30; %densité critique (en veh/m)
Vhat=22; %vitesse critique (en m/s)
beta=(V-Vhat)/Rhat; %paramètre de pente de la fonction vitesse en fluide
qstar=1; %paramètre

A=V+qstar/R-I;
B=beta-I/R;
rhocrit=1/(2*B)*(A-sqrt(A^2-4*qstar*B));%densité critique
rcrit=1/rhocrit;

if mu==0
    L=Colombo(inf,I);
else %mu >0
    v1=sqrt(beta/mu);

    g1=@(r) V-beta/r-mu*r;
    g2=@(r) (I+qstar*r)*(1-1/(R*r))-mu*r;

    if rcrit <= v1
        if mu<=qstar
            L=max(g1(v1),g2(rcrit));
        else %mu>qstar
            v2=sqrt(I/(R*(mu-qstar)));
            if v2>rcrit
                L=max(g1(v1),g2(rcrit));
            elseif v2<=rcrit && v2>=rmin
                L=max(g1(v1),g2(v2));
            elseif v2<rmin
                L=max(g1(v1),g2(rmin));
            end
        end
    else %rcrit > v1
        if mu<=qstar
            L=max(g1(rcrit),g2(rcrit));
        else %mu>qstar
            v2=sqrt(I/(R*(mu-qstar)));
            if v2>rcrit
                L=max(g1(rcrit),g2(rcrit));
            elseif v2<=rcrit && v2>=rmin
                L=max(g1(rcrit),g2(v2));
            elseif v2<rmin
                L=max(g1(rcrit),g2(rmin));
            end
        end
    end
end
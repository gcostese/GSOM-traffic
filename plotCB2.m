function plotCB2(t0,tmax,Deltat,Xn0,vn0,Rn0)
%fonction qui permet de tracer les conditions au bord

figure
h = plot((t0:Deltat:tmax),Xn0,'b','Linestyle','-','LineWidth',2);
xlabel('Time t (s)','Fontsize',16)
ylabel('Position X_0 (m)','Fontsize',16)
title('Position X(N_0,t)','Fontsize',18)
saveas(h,'initialposition.eps','epsc')

V=@(t) vn0(fix((t-t0)/Deltat)+1);
figure
h3 = ezplot(V,[t0 tmax]);
set(h3, 'Color','r','Linestyle','-','LineWidth',2);
xlabel('Time t (s)','Fontsize',16)
ylabel('Speed v_0 (m.s^{-1})','Fontsize',16)
title('Speed v(N_0,t)','Fontsize',18)
saveas(h3,'vbound.eps','epsc')

R=@(t) Rn0(fix((t-t0)/Deltat)+1);
figure
h4 = ezplot(R,[t0 tmax]);
set(h4, 'Color','r','Linestyle','-','LineWidth',2);
xlabel('Time t (s)','Fontsize',16)
ylabel('Spacing r_0 (m.s^{-1})','Fontsize',16)
title('Spacing r(N_0,t)','Fontsize',18)
saveas(h4,'rbound.eps','epsc')
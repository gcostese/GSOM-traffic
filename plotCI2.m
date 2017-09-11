function plotCI2(N0,Nmax,DeltaN,Rt0,I,Xt0)
%fonction qui permet de tracer les conditions initiales r_0(N) et I_0(N)

R=@(n) Rt0(fix((n-N0)/DeltaN)+1);
figure
h = ezplot(R,[N0 Nmax]);
set(h, 'Color','r','Linestyle','-','LineWidth',2);
xlabel('Label N','Fontsize',16)
ylabel('Spacing r_0 (m)','Fontsize',16)
title('Initial conditions r(N,t_0)','Fontsize',18)
saveas(h,'rini.eps','epsc')

figure
h2 = ezplot(I,[N0 Nmax]);
set(h2, 'Color','b','Linestyle','-','LineWidth',2);
xlabel('Label N','Fontsize',16)
ylabel('Driver attribute I_0','Fontsize',16)
title('Initial conditions I(N,t_0)','Fontsize',18)
saveas(h2,'Iini.eps','epsc')

figure
h3 = plot((N0:DeltaN:Nmax),Xt0,'b','Linestyle','-','LineWidth',2);
xlabel('Label N','Fontsize',16)
ylabel('Position X (m)','Fontsize',16)
title('Initial positions X(N,t_0)','Fontsize',18)
saveas(h3,'Xini.eps','epsc')
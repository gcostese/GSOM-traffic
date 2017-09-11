function [g,S]=plotSpeed(X,Deltat,t0,tmax,N0,Nmax)
%Function that allows to draw the values of speed over the computational
%domain [N0 Nmax]*[t0 tmax]

[n,m]=size(X);
S=zeros(n-1,m);

for i=2:n
    for j=1:m
        S(i,j)=(X(i,j)-X(i-1,j))/Deltat;
    end
end

%Parameter
vmax=25; %in m/s

g=figure;
image('xdata',[N0 Nmax],'ydata',[t0 tmax],'cdata',S)
axis([N0 Nmax t0 tmax])
colormap(jet(vmax))
colorbar
title('Speeds (m/s)','FontSize',18)
xlabel('Label','FontSize',16)
ylabel('Time (s)','FontSize',16)
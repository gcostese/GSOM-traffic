function S=Speed(v1,v2,t0,T,Deltat)
%Function that allows to compute the trajectory of a simulated vehicle
%between times t0 and T, according to variations of time length Deltat
%First, its acceleration is computed thanks to a centered white noise with
%a variance of 2 m/s^2, and then the speed is computed as the integral of
%the acceleration over time
%The profile of mean speed is piecewise constant, with a first value aroud
%v1 (between t0 and middle of time span) and a second value v2 (from the
%middle of time span to T)

% close all

lambda=1;
sigma=3; %en m/s

n=floor((T-t0)/Deltat)+1;
S=zeros(1,n);

A=sqrt(sigma*sqrt(Deltat)).*randn(1,n);
% %Vérification graphique de la répartition de A
% figure;
% hist(A,100)
% disp(abs(var(A)-sigma*sqrt(Deltat)))
% f=@(t) A(floor((t-t0)/Deltat)+1);
% figure;
% fplot(f,[t0 T]);
% title('White noise')


%Initialization
w=0;
ind=1;
for t=t0:Deltat:(T+t0)/2
    w=(1-lambda*Deltat)*w+A(ind);
    S(ind)=v1+w;
    ind=ind+1;
end
for t=(T+t0)/2+Deltat:Deltat:T
    w=(1-lambda*Deltat)*w+A(ind);
    S(ind)=v2+w;
    ind=ind+1;
end

% %Vérification graphique
% g=@(t) S(floor((t-t0)/Deltat)+1);
% h=@(t) v1.*(t<=t0+(T-t0)/2)+v2.*(t>t0+(T-t0)/2);
% figure;
% hold on
% fplot(g,[t0 T]);
% fplot(h,[t0 T],'r');
% title('Speed profile')
% hold off
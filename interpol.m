function v=interpol(x,a,b,c,va,vb,vc)
%fonction d'interpolation dans un triangle (a,b,c) où a, b et c sont les
%vecteurs colonnes des coordonnées dans le plan
%x est le vecteur colonnes des coordonnées du point où on cherche à
%interpoler la valeur v
%la valeur v en x est déterminée en fonction des valeurs à chaque sommet du
%triangle

n=5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if va==inf || vb==inf || vc==inf
%     disp('valeur infinie dans l''interpolation!');
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A=[a b c; 1 1 1];
B=[x; 1];

%calcul des coordonnées barycentriques du point x
l=A\B; %equivalent à inv(A)*B

%calcul de la valeur v
if round(l(1)*10^n)>=0 && round(l(2)*10^n)>=0 && round(l(3)*10^n)>=0
    v=l(1)*va+l(2)*vb+l(3)*vc;
else
    v=inf;
end

%sécurité pour éviter que le min attrape une valeur infinie (négativement)
if v==-inf
    v=inf;
end
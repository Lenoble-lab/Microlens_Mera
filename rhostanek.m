function r = rhostanek(R,z,th)

% densite massique du bulbe : modele de Stanek
% valeurs selon S. Calchi Novati et al., Astron.Astrophys. 480 (2008) 723, 0711.3758.

global Rcoro

%Rcoro est le rayon de corotation
%la densité est tronquée pour r>Rcoro
%Dans l'article, Rcoro = 3.5 kpc

rho0=1.0;       % valeur test

i0 = find(R>Rcoro);
i1 = find(R<=Rcoro);

%%-----------------------%%
%%Modèle de Calchi Novatti 2008, Stanek et al 1997
%%-----------------------%%
phi=23.8; %angle de la barre avec l'axe soleil-GC
alpha=pi*phi/180;
beta=0;
x0=890;
y0=x0*4.3/10;
z0=x0*2.8/10;


% calcul des coordonnees dans le repere principal du bulbe

ca=cos(alpha);
sa=sin(alpha);
cb=cos(beta);
sb=sin(beta);

cth=cos(th);
sth=sin(th);



x=R.*cth*ca+R.*sth*sa;
y=-R.*cth*sa*cb+R.*sth*ca*cb+z*sb;
z=R.*cth*sa*sb-R.*sth*ca*sb+z*cb;



a=sqrt((x/x0).^2+(y/y0).^2+(z/z0).^2);

r(i0) = zeros(size(i0))
r(i1) = rho0*exp(-a);
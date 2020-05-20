function r = rhostanek(R,z,th)

% densite massique du bulbe : modele de Stanek
% valeurs selon S. Calchi Novati et al., Astron.Astrophys. 480 (2008) 723, 0711.3758.

global Rcoro

%Rcoro est le rayon de corotation
%la densité est tronquée pour r>Rcoro
%Dans l'article, Rcoro = 3.5 kpc



i0 = find(R>Rcoro);
i1 = find(R<=Rcoro);

%%-----------------------%%
%% Modèle de Calchi Novatti 2008, Stanek et al 1997
%%-----------------------%%
phi=28; %angle de la barre avec l'axe soleil-GC
alpha=pi*phi/180;
beta=0;
x0=890;
y0=x0.*4.3/10;
z0=x0.*2.8/10;


M_b = 1.8e10;
rho0 = M_b/(x0*y0*z0*8*pi);

% calcul des coordonnees dans le repere principal du bulbe

ca=cos(alpha);
sa=sin(alpha);
cb=cos(beta);
sb=sin(beta);

cth=cos(th);
sth=sin(th);


x=R.*cth*ca+R.*sth*sa;
y=-R.*cth*sa*cb+R.*sth*ca*cb+z.*sb;
z=R.*cth*sa*sb-R.*sth*ca*sb+z.*cb;

a = sqrt((x./x0).^2+(y./y0).^2+(z./z0).^2);

r = rho0.*exp(-a);

%Iocco (avec une coupure exponentielle
% r(i0) = r(i0).*exp(-(R(i0)-Rcoro).^2/(2*500^2));

%Calchi Novatti (arrêt brusque)
% r(i0) = zeros(size(i0));

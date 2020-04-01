
% pour l'integration


% Densite du buble, modele de Zhao

% Parametres : coordonnees cylindriques galactiques en pc R et z, et theta
% retourne la densite de masse en Msol/pc^3
%
% alpha est l'angle entre l'axe soleil-centre et l'axe de la barre
% beta est l'angle entre l'axe z et le nouvel axe z de l'ellipsoide
% (ie angle de rotation pour passer des axes de la galaxie y,z aux axes
% principaux de la barre
%

function rh = rhozhao(R,z,th)

%--------------------------
% donnees du modele
% a modifier eventuellement
%--------------------------

%rho0=1.1968;         % essai

%rho0=1.8469;         % masse du bulbe = 1.5e10 
%rho0=2.71;         % masse du bulbe = 2.2e10 
rho0=1.2313;       % masse du bulbe = 1.0e10
%rho0=2.4625;         % masse du bulbe = 2.e10
%rho0=4;

alpha=pi*20/180;
beta=0;
x0=1490;
y0=580;
z0=400;
qa=0.6;

%------------------------
%intermediaires de calcul
%------------------------


alpha=alpha;
ca=cos(alpha);
sa=sin(alpha);
cb=cos(beta);
sb=sin(beta);

cth=cos(th);
sth=sin(th);

X=R.*cth*ca+R.*sth*sa;
Y=-R.*cth*sa*cb+R.*sth*ca*cb+z*sb;
Z=R.*cth*sa*sb-R.*sth*ca*sb+z*cb;

sb2=sqrt(((X/x0).^2+(Y/y0).^2).^2+(Z/z0).^4);
sa=sqrt((qa*qa*(X.^2+Y.^2)+Z.^2)/(z0^2));

%-----------------
% densite de masse
%-----------------

%rh=rho0*(exp(-sb2/2));   % modele de Dwek (ref du modele G2)

rh=rho0*(exp(-sb2/2)+sa.^(-1.85).*exp(-sa)).*R;

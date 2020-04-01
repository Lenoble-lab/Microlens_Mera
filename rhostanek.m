function r = rhostanek(R,z,th)

% densite massique du bulbe : modele de Stanek



%rho0=;         % masse du bulbe = 2.2e10 
rho0=1.0;       % masse du bulbe = 1.0e10


phi=25; %angle phi compris entre 20 et 30 degre
alpha=pi*phi/180;
beta=0;
x0=900;
y0=x0*4/10;
z0=x0*3/10;

% calcul des coordonnees dans le repere principal du bulbe

alpha=alpha;
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

r = rho0*exp(-a);
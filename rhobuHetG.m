% Densite du buble, modele utilisé par Han&Gould (1995 et 2003), également utilisé dans Calchi Novatti 2008
    
% Parametres : coordonnees cylindriques galactiques en pc R et z, et theta
% retourne la densite de masse en Msol/pc^3
%
% alpha est l'angle entre l'axe soleil-centre et l'axe de la barre
% beta est l'angle entre l'axe z et le nouvel axe z de l'ellipsoide
% (ie angle de rotation pour passer des axes de la galaxie y,z aux axes
% principaux de la barre
%

function res = rhobuHetG(R,z,th)

    global Ro 

    s = (R.^4 + (z./0.610).^4).^(1/4);

    i0 = find(s<938);
    i1 = find(s>=938 & R<700);
    i2 = find(R>=700);


    res(i0) = 1.04e6.*(s(i0)./0.482).^(-1.85);
    res(i1) = 3.53 * besseli(0, s(i1)/667);

    rho0 = 3.66e-2; %en M_sol/pc^3
    
    
    
    %%-----------------------%%
    %%Modèle en dehors de la zone centrale, à l'origine le modèle G2 de Dwek 1995
    %%-----------------------%%
    alpha=pi*20/180;
    beta=0;
    x0=1580; %en parsec
    y0=620;
    z0=430;

    
    %------------------------
    %intermediaires de calcul
    %------------------------
    
    
    ca=cos(alpha);
    sa=sin(alpha);
    cb=cos(beta);
    sb=sin(beta);
    
    cth=cos(th);
    sth=sin(th);
    
    X=R.*cth*ca+R.*sth*sa;
    Y=-R.*cth*sa*cb+R.*sth*ca*cb+z*sb;
    Z=R.*cth*sa*sb-R.*sth*ca*sb+z*cb;
    
    sb2=sqrt(((X(i2)/x0).^2+(Y(i2)/y0).^2).^2+(Z(i2)/z0).^4);
    %sa=sqrt((qa*qa*(X.^2+Y.^2)+Z.^2)/(z0^2));
    
    %-----------------
    % densite de masse
    %-----------------
    
    res(i2)=rho0*(exp(-sb2/2));  
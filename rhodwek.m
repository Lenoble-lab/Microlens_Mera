% Densite du buble, modele de Dwek (ref du modele G2)
    
%
% Parametres : coordonnees cylindriques galactiques en pc R et z, et theta
% retourne la densite de masse en Msol/pc^3
%
% alpha est l'angle entre l'axe soleil-centre et l'axe de la barre
% beta est l'angle entre l'axe z et le nouvel axe z de l'ellipsoide
% (ie angle de rotation pour passer des axes de la galaxie y,z aux axes
% principaux de la barre
%

function rh = rhodwek(R,z,th)

    global Ro Rcoro
    %--------------------------
    % donnees du modele
    % a modifier eventuellement
    %--------------------------
    
%     rho0 = 1;
   %-------------------------
   %Donnee de calchi novatti
   %------------------------
%     rho0 = 2.4 ; %Masse du bulbe 1.410^10M_sol (Calchi Novatti)

   %-------------------------
   %Donnee de Iocco
   %------------------------
%     rho0 = 2.4 * 0.73; %Iocco, modèle 2

    rho0 = 2.4 *0.94;

    i0 = find(R>Rcoro);
    i1 = find(R<=Rcoro);
    
    %%-----------------------%%
    %%Modèle de Iocco (G2), 2018
    %%-----------------------%%
    alpha=pi*24.9/180;
    beta=0;
    x0=1239 * Ro/8000;; %en parsec
    y0=609 * Ro/8000;;
    z0=438 * Ro/8000;;
    qa=0.6 * Ro/8000;;
    
    
    %%-----------------------%%
    %%Modèle de Han&Gould, Calchi Novatti 2008, à l'origine le modèle G2 de Dwek 1995
    %%-----------------------%%
    % alpha=pi*20/180;
    % beta=0;
    % x0=1580; %en parsec
    % y0=620;
    % z0=430;

    %%-----------------------%%
    %%Valeures présentes dans le code à l'origine
    %%-----------------------%%
    % alpha=pi*20/180;
    % beta=0;
    % x0=1490; %en parsec
    % y0=580;
    % z0=400;

    
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
    
    sb2=sqrt(((X/x0).^2+(Y/y0).^2).^2+(Z/z0).^4);
    %sa=sqrt((qa*qa*(X.^2+Y.^2)+Z.^2)/(z0^2));
    
    %-----------------
    % densite de masse
    %-----------------
    
    rh=rho0*(exp(-sb2/2));  

    %coupure exponentielle (Iocco)
    rh(i0) = rh(i0).*exp(-(R(i0)-Rcoro).^2/(2*500^2));

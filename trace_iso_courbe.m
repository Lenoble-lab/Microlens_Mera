%tracé des isocourbe de tau selon l'angle (l,b) en direction du bulbe


%----------------------------------
% Constantes physiques (unites SI)
%----------------------------------
G=6.672e-11;	pc=3.08567802e16;
kpc=pc*1e3;  	Msol=1.989e30;
c=299792458;	GMsol=1.32712497e20;


%----------------------------------------------------------------------
% Param�tres de la fonction de distribution de la distance de la source
%----------------------------------------------------------------------

global dsup dinf 

dsup = 10000.;
dinf = 800.; %distance en parsec

%-------
% soleil 
%-------
global Ro elev
Ro = 8000;	   % distance Sun-GC en pc
elev = 0 ;      % elevation au dessus du plan du disque
Lkpc=Ro/1000;  % distance Sun-GC en kpc


%----------------------------
%rayon de corotation (en pc)
%----------------------------

global Rcoro
Rcoro = 3500;



uT = 1;	% donnée ogle	   % Seuil de d�tection en param�tre d'impact
AT = 3/sqrt(5);    % Seuil de d�tection en amplification


L = (-1:1e-5:1).*10;
B = (-1:1e-5:1).*10;

tau_table = tau(L, B)


figure(1)
hold on
contour(L, B, tau_table, 'ShowText','on')
hold off

%-------------------------------
%calcul de la profondeur optique, let b en degré
%-------------------------------

function res = tau(l1, b1)
    %---------------------------------------       
    % param�tres du calcul de microlentilles
    %---------------------------------------
    global l b

    % definition de la fenetre de Baade :  l = 1 et b = -4
    l = l1 *pi/180;    % direction d'observation en radian
    b = b1 *pi/180;

    %-------------------------------------------------
    % Calcul pr�alable du cosinus pour aller plus vite
    %-------------------------------------------------
    global sinb cosb  cosbl sinl cosl
    sinb = abs(sin(b));		cosb = cos(b);		cosl = cos(l);
    cosbl=cos(b)*cos(l);		sinl = sin(l);


    normnu=quadl('nsource',dinf,dsup);
    Ctau = 4*pi*GMsol*uT*uT/c/c/pc;
    tau  = Ctau*dblquad('dtau',0,1,dinf,dsup,1e-2,'quadl') /normnu;
    res=real(tau);
end


% Fonction a integrer pour le calcul de profondeur optique 

function res = dtau(x,L)
    global Ro cosb cosbl sinb bet
    
    res = rho_tot(x.*L).*x.*(1-x).*nsource(L).*L.^2;
end

%% Nombre de source 
%%Prend en compte la fonction de luminosité
function res = nsource (x)

    global dsup dinf 
    
    res = zeros(size(x));
    
    bet=-1;
    
    
    i1 = find(x>=dinf & x<=dsup);
    if (length(i1)>=1)
      res (i1) = rho_tot(x(i1)).*x(i1).^(2.*bet+2);
    end
end

%%densité totale bulbe+disque
function res = rho_tot(x)

    [R, z, th] = toGC(x);
    
    res=zeros(size(R));
    
    %---------------
    % source : bulbe
    %---------------
    
    res = res + rhostanek(R,z,th);
    
    %----------------------
    % source : disque 
    %----------------------
    
    res = res + rhodHetG(R,z,th);
    
end
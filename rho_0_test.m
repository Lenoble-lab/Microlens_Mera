%% but du prgm : déterminer la densité rho_0 du bulbe à partir de différents modèles et des observations
%% de profondeur optique en direction du bulbe

%% Principe : on part de rho_0 = 1 et ensuite comme \tau en dépend, on le trouve pour que la valeur calculée correspondent à la valeur observée


%----------------------------------
%% Constantes physiques (unites SI)
%----------------------------------
G=6.672e-11;	pc=3.08567802e16;
kpc=pc*1e3;  	Msol=1.989e30;
c=299792458;	GMsol=1.32712497e20;


%----------------------------------------------------------------------
%% Param�tres de la fonction de distribution de la distance de la source
%----------------------------------------------------------------------

global dsup dinf 

dsup = 30000.;
dinf = 0.; %distance en parsec

%-------
% soleil 
%-------
global Ro elev
Ro = 8000;	   % distance Sun-GC en pc
elev = 0 ;      % elevation au dessus du plan du disque
Lkpc=Ro/1000;  % distance Sun-GC en kpc


%----------------------------
%%rayon de corotation (en pc)
%----------------------------

global Rcoro
Rcoro = 3500;


%---------------------------------------       
%% param�tres du calcul de microlentilles
%---------------------------------------
global l b

% definition de la fenetre de Baade :  l = 1 et b = -4
l = 1.5 *pi/180;    % direction d'observation en radian
b = -2.68 *pi/180;

uT = 1;	% donnée ogle	   % Seuil de d�tection en param�tre d'impact
AT = 3/sqrt(5);    % Seuil de d�tection en amplification



%-------------------------------------------------
% Calcul pr�alable du cosinus pour aller plus vite
%-------------------------------------------------
global sinb cosb  cosbl sinl cosl
sinb = abs(sin(b));		cosb = cos(b);		cosl = cos(l);
cosbl=cos(b)*cos(l);	sinl = sin(l);

%-------------------------------
%% calcul de la profondeur optique
%-------------------------------

x = (0:1e-3:1).*(dsup-dinf)+dinf;

tau = zeros(size(x));
y_2 = zeros(size(x));

% for k = 1:numel(x)
%     y(k) = integral(@nsource_1,0,x(k));
%     y_2(k) = quadl(@nsource_1,0,x(k));
% end


dsup = 20000;
dinf = 0;

normnu=integral(@nsource_1,dinf,dsup);
Ctau = 4*pi*GMsol*uT*uT/c/c/pc;
tau  = Ctau*integral2(@dtau_1,0,1,dinf,dsup) /normnu;
taucx=tau;
tau=real(taucx);
disp(['nsource = ' num2str(normnu*1e-11)]);
disp(['tau = ' num2str(tau)]);


% Fonction a integrer pour le calcul de profondeur optique 

function res = dtau_1(x,L)
    
    res = rho_tot(x.*L).*x.*(1-x).*nsource_1(L).*L.^2;
end

%% Nombre de source 
%%Prend en compte la fonction de luminosité
function res = nsource_1 (x)

    global dsup dinf 
    
    res = zeros(size(x));
    
    bet=0;
    
    
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
    res = res + rhostanek(R, z, th);
%     res = res + rhodwek(R, z, th);
    
    %----------------------
    % source : disque 
    %----------------------
    
    res = res + rhodHetG(R,z,th);
    
end
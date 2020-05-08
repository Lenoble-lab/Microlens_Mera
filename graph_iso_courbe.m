%tracé des isocourbe de tau selon l'angle (l,b) en direction du bulbe

fichier='graph_iso.mat';

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
%rayon de corotation (en pc)
%----------------------------

global Rcoro
Rcoro = 3500;



uT = 1;	% donnée ogle	   % Seuil de d�tection en param�tre d'impact
AT = 3/sqrt(5);    % Seuil de d�tection en amplification


L = (0:1e-2:1).*20-10;
B = (0:1e-2:1).*20-10;


% L = [1.5, 1.16, -4.5, -1.5, 1.5, 4.5, 306.56, 331.09, 18.51, 26.6];
% B = [-2.68, -2.75, 2.4, 2.42, 2.22, 2.53, -1.46, -2.42, -2.09, -2.15];


% L = [1.5, 6, 8.5, 4];
% B = [-2.68, 0, 0, 4];


Ctau = 4*pi*GMsol*uT*uT/c/c/pc;

tau_table = zeros(size(L,2), size(B,2));
% tau_table = zeros(size(L,2));

for k = 1:numel(L)
        
    if ( rem(k,1) == 0 )
    
disp(['compteur = ' num2str(k)]);

    end
%     tau_table(k) = Ctau * tau(L(k), B(k));
%     disp(['l = ' num2str(L(k)) '  b = ' num2str(B(k)) '  tau = ' num2str(tau_table(k))])
    
    for i = 1:numel(B)
        tau_table(k, i) = Ctau * tau(L(i), B(k));
        
    end
end

disp([' '])

save('graph_iso.mat', 'tau_table', 'B', 'L')

%% test densité
R = (0:1e-3:1).*4000;
th = (0:1e-3:1).*2*pi;

[R, th] = meshgrid(R, th);
rho_bu = rhostanek(R, zeros(size(R)), th);
rho_disk = rhodHetG(R, zeros(size(R)), th);
[X,Y] = pol2cart(th, R);

rho = rho_disk  ;

figure(1)
pcolor(X, Y, rho); 
shading interp;
hold on
colorbar;


contour(X, Y, rho, '--black', 'ShowText', 'on');

%% analyse the results


clear
tau_load = load('graph_iso.mat');

figure(1)
hold on
% contour(tau_load.B, tau_load.L, transpose(tau_load.tau_table), 'ShowText', 'on')
contour(tau_load.B, tau_load.L, transpose(tau_load.tau_table), [0.3e-6, 1e-6, 1.5e-6, 2.3e-6, 3.5e-6], 'ShowText', 'on')
% contour(tau_load.B, tau_load.L, transpose(tau_load.tau_table), [3.5e-6, 4e-6, 1.5e-6, 2.17e-6, 3e-6], 'ShowText', 'on')
% contour(tau_load.B, tau_load.L, fliplr(transpose(tau_load.tau_table)),[3.5e-6, 6e-6, 1.5e-6, 2.17e-6, 3e-6], 'ShowText', 'on')


%-------------------------------
%calcul de la profondeur optique, let b en degré
%-------------------------------

function res = tau(l1, b1)
    %---------------------------------------       
    % param�tres du calcul de microlentilles
    %---------------------------------------
    global l b

    l = l1 *pi/180;    % direction d'observation en radian
    b = b1 *pi/180;

    %-------------------------------------------------
    % Calcul pr�alable du cosinus pour aller plus vite
    %-------------------------------------------------
    global sinb cosb  cosbl sinl cosl
    sinb = abs(sin(b));		cosb = cos(b);		cosl = cos(l);
    cosbl=cos(b).*cos(l);		sinl = sin(l);
    
    global dsup dinf

    normnu=integral(@nsource,dinf,dsup);
    taux  = integral2(@dtau,0,1,dinf,dsup, 'Method', 'iterated') /normnu;
    res=real(taux);
end


% Fonction a integrer pour le calcul de profondeur optique 




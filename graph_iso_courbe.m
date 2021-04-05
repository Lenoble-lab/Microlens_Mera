%tracé des isocourbe de tau selon l'angle (l,b) en direction du bulbe
clear
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

global vlimit

vlimit = 1000e3;

%----------------------
% Nombre de simulations
%----------------------

n = 100*1e5;
% n = 5000;
nbsimul=1; %a augmenter pour meilleure stat

%----------------------------------------------------------------------
% Param�tres de la fonction de distribution de la distance de la source
%----------------------------------------------------------------------

global dsup dinf 

dsup = 20000.;
dinf = 800.;
%distance en parsec

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

%-----------------------------------
% Param�tres de la fonction de masse
%-----------------------------------
global minf msup

minf=0.01;
msup=100;  

global Vinf Vsup 

Vinf = 22; % magnitudes limites des étoiles observées.
Vsup = 16;

uT = 1;	% donnée ogle	   % Seuil de d�tection en param�tre d'impact
AT = 3/sqrt(5);    % Seuil de d�tection en amplification


L_table = ((1e-5:0.05:1).*540-270)*pi/180;
B_table = ((1e-5:0.2:1).*14-7)*pi/180;

Ctau = 4*pi*GMsol*uT*uT/c/c/pc;

global l b
        
tau_table = zeros(size(L_table,2), size(B_table,2));
gamma_table = zeros(size(L_table,2), size(B_table,2));

for l_indice = 1:numel(L_table)   
    for b_indice = 1:numel(B_table)
    
    %On definti les angles : 
    l = L_table(l_indice); 
    b = B_table(b_indice);
    %calcul avec les angles indiqué
    main
    
    %Valeur de tau
    tau_table(l_indice, b_indice) = tau;
   
    %----------------------------------------
% recuperation des evenements selectionnes
%----------------------------------------

load evenements.txt
x=evenements(:,1) ; x=x';
ds=evenements(:,2) ; ds=ds';
v=evenements(:,3) ; v=v';
m=evenements(:,4) ; m=m';
te=evenements(:,5) ; te=te';
    
    %Valeur de gamma
    gam1=4*sqrt(GMsol)/c*uT/sqrt(pc*pc*pc)*length(te)/(n*nbsimul)*86400*365.25*1e6/mmean;
    gamma_table(l_indice, b_indice) = gam1*Gammax;
        
    end
end

disp([' '])

save('graph_iso_model.mat', 'tau_table','gamma_table', 'B_table', 'L_table')


%% analyse the results
clear;close all;
tau_load = load('graph_iso_model.mat');

% Convertion radian->degré
tau_load.B_table = tau_load.B_table * 180/pi;
tau_load.L_table = tau_load.L_table * 180/pi;

[L_table, B_table] = meshgrid(tau_load.B_table, tau_load.L_table);



figure(2)
pcolor(L_table, B_table, log10(tau_load.tau_table)); 
% shading interp;
hold on
colorbar;
% contour(x_y, y_x, log10(dens_surf), 'black', 'ShowText', 'on');
xlabel('longitude galactique (en degrès)')
ylabel('latitude galactique (en degrés)')

figure(1)
hold on
% open_table7_1d
% contour(tau_load.B, tau_load.L, transpose(tau_load.tau_table), 'ShowText', 'on')
% contour(tau_load.B, tau_load.L, tau_load.tau_table.*1e6, [0.3, 1, 1.5, 2.17, 3], 'ShowText', 'on')
contour(tau_load.L_table, tau_load.B_table, transpose(tau_load.tau_table).*1e6, [0.03, 0.06, 0.09, 0.3, 1, 1.5, 2.17, 3], 'ShowText', 'on')

% contour(tau_load.B, tau_load.L, transpose(tau_load.tau_table), [3.5e-6, 4e-6, 1.5e-6, 2.17e-6, 3e-6], 'ShowText', 'on')
% contour(tau_load.B, tau_load.L, fliplr(transpose(tau_load.tau_table)),[3.5e-6, 6e-6, 1.5e-6, 2.17e-6, 3e-6], 'ShowText', 'on')
axis('equal')
grid on;
grid minor
xlabel('longitude galactique (en degrès)')
ylabel('latitude galactique (en degrés)')


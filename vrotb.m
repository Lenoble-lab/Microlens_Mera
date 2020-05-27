    % Courbe de rotation pour les etoiles du bulbe
%
% Parametre : les coordonnees galactiques
% retourne la vitesse vphi en m/s

function v=vlrot(R,z,th)

global Rcoro

%v=ones(size(R)).*100e3; 

%--------------------------------------
% rotation rigide, omega = 39 kms/s/kpc Portail et al, 2017 
%--------------------------------------

v=39.1e3.*R./1000;   

%--------------------------------------
% rotation rigide, omega = 60 kms/s/kpc 
%--------------------------------------

% v=57.1e3.*R./1000;    % tire de englemier et gerhard


%-------------------------
% rayon de corotation a Rc
%-------------------------

%v = 200e3.*R./Rcoro;           % Rcoro doit etre en pc
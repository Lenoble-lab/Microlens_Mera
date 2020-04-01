% Densite du disque epais
%
% Parametres : coordonnees cylindriques galactiques en pc R et z, et theta
% retourne la densite de masse en Msol/pc^3

function res = rhotd(R,z,th)
global Ro 

% ------------
% modele de KP
% ------------

%res = rhoKP(R,z,th);

%-------------------
% disque double expo
%-------------------

htd = 760;	        % Echelle de hauteur en pc
Ltd = 3000; 	        % Echelle de longueur en pc
rhtd = 0.05/10;        % densite au voisinage solaire, en Msol/pc^3

res = rhtd.*exp(-R./Ltd-abs(z)./htd+Ro./Ltd);

%--------------------
% pas de disque epais
%--------------------

%res = zeros(size(R));





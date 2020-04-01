% Courbe de rotation pour les etoiles du disque
%
% Parametre : les coordonnees galactiques
% retourne la vitesse vphi en m/s

function v=vlrot(R,z,th)
v=ones(size(R)).*200e3;  	

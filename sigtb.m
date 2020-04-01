% Dispersion de vitesse des etoiles du bulbe dans la direction theta
%
% Parametres : coordonnees Galactiques en pc
% Sortie: dispersion de vitesse en m/s

function res = siglp(R,z,t)
res = ones(size(R)).*110e3*sqrt(2);
%res = ones(size(R)).*140e3;
% Dispersion de vitesse des etoiles du bulbe dans la direction theta
%
% Parametres : coordonnees Galactiques en pc
% Sortie: dispersion de vitesse en m/s

function res = sigtb(R,z,t)
res = ones(size(R)).*110e3;
%res = ones(size(R)).*140e3;
% res = zeros(size(R));
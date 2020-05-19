% Dispersion de vitesse des etoiles du bulbe dans la direction z
%
% Parametres : coordonnees Galactiques en pc
% Sortie: dispersion de vitesse en m/s

function res = sigzb(R,z,t)
res = ones(size(R)).*110e3;
% res = ones(size(R)).*120e3;

% res = zeros(size(R));
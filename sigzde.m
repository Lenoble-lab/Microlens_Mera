% Dispersion de vitesse des etoiles du disque dans la direction z
%
% Parametres : coordonnees Galactiques en pc
% Sortie: dispersion de vitesse en m/s

function res = sigtde(R,z,t)
    %res = ones(size(R)).*20e3;
    
    %-----------
    %Pasetto et al
    %-----------
    
    res = ones(size(R)).*35.1;
    
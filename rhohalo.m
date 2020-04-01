% Densite de masse du halo sombre
% Parametres : coordonnees cylindriques galactiques en pc R et z et theta
% retourne la densite de masse en Msol/pc^3

function res = rhoh(R,z,th)
global Ro 


%-----------------
% modele isotherme
%-----------------

rh0 = 7.9e-3;	        % normalisation du halo au voisinage solaire
a = 5000;		% rayon de coeur du halo en pc
q = 1;	                % aplatissement du halo

%res = rh0.*(Ro*Ro+a*a)./(R.*R+z.*z./q^2+a*a);

%------------
% pas de halo
%------------

res = zeros(size(R));

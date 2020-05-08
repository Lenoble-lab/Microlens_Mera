% Densite du buble
%
% Parametres : coordonnees cylindriques galactiques en pc R et z, et theta
% retourne la densite de masse en Msol/pc^3


function res = rhobulbe(R,z,th)

%---------------
% modele de Kent
%---------------

%res = rhobk(R,z,th);

%-----------------
% modele de Dweck
%-----------------

% res = rhodwek(R,z,th);

%-----------------
% modele de Stanek
%-----------------

res = rhostanek(R,z,th);

%---------------
% modele de Zhao
%---------------

% res = rhozhao(R,z,th);

%-----------------------
% modele de Zhao tronque
%-----------------------

%res = rhozhaotrunc(R,z,th);

%---------------
% modele de Bale
%---------------

%res = rhoBGSb(R,z,th);


%-------------
% pas de bulbe
%-------------


%res = zeros(size(R));
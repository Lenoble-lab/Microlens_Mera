% Densite de masse sous forme d'etoiles dans le disque
%
% Parametre : coordonnees cylindriques galactiques en pc R et z, et theta
% Retourne la densite en Msol/pc^3
%
%

function res = rhodm(R,z,th)

global Ro

% ------------
% modele de KP
% ------------

%res = rhoKP(R,z,th);

%-------------------
% disque double expo
%-------------------

hd=250;
Ld=3000;
rhsol=0.05;

% res = rhsol*exp(-R./Ld-abs(z)./hd+Ro./Ld);
res = rhodHetG(R,z,th);

%--------------------
% pas de disque mince
%--------------------

%res = zeros(size(R));

%--------------
% disque de BGS
%--------------


%res = rhoBGSd(R,z,th);

% Densite de lentilles
%
% Parametre : coordonnees cylindriques galactiques en pc R et z, et theta
% Retourne la densite en Msol/pc^3
%

function res = rholens (x)

global mmeandm mmeanbu mmeande mmeanh

[R, z, th] = toGC(x);

res=zeros(size(R));

%---------------
% lentille : bulbe
%---------------

res = res + rhobulbe(R,z,th);

%----------------------
% lentille : disque H et G
%----------------------

res = res + rhodHetG(R, z, th);

%----------------------
% lentille : disque mince
%----------------------

% res = res + rhodm(R,z,th);

%----------------------
% lentille : disque epais
%----------------------

% res = res + rhode(R,z,th);

%--------------
% lentille : halo
%--------------

%res = res + rhohalo(R,z,th);


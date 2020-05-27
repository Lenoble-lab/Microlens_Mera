% Dispersion de vitesse des etoiles du disque dans la direction z
% orthogonale au plan
%
% Parametres : coordonnees Galactiques en pc
% Sortie: dispersion de vitesse en m/s

function res = sigzdm(R,z,t)

res = ones(size(R)) * 20e3;

global Rcoro Ro


% donnée de Pasetto, S. et al 
%"Thin disk kinematics from RAVE and the solar motion", Astronomy and Astrophysics 547 (2012), pp. A71.

% res = ones(size(R)) * 16.3e3;
% res = zeros(size(R));

% i0 = find( R <= Rcoro );   
% i1 = find( R > Rcoro & R <= (Rcoro+Ro)/2 );
% i2 = find( R > (Rcoro+Ro)/2 );
% 
% p=(20e3-110e3)*sqrt(2)/(Ro-Rcoro);
% res(i1)=p*(R(i1)-Rcoro)+110e3*sqrt(2);
% 
% res(i2)=ones(size(i2)).*20e3*sqrt(2);
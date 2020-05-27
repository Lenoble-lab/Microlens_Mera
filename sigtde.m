% Dispersion de vitesse des etoiles du disque dans la direction theta
%
% Parametres : coordonnees Galactiques en pc
% Sortie: dispersion de vitesse en m/s

function res = sigtde(R,z,t)
%res = ones(size(R)).*20e3;

%-----------
%Pasetto et al
%-----------

res = ones(size(R)).*46.1;

%-----------
%Robin et al
%-----------

%res = ones(size(R)).*20e3*sqrt(2);
%res = ones(size(R)).*20e3;


% global Rcoro Ro


% i0 = find( R <= Rcoro );   
% i1 = find( R > Rcoro & R <= (Rcoro+Ro)/2 );
% i2 = find( R > (Rcoro+Ro)/2 );

% p=(20e3-110e3)*sqrt(2)/(Ro-Rcoro);
% res(i1)=p*(R(i1)-Rcoro)+110e3*sqrt(2);

% res(i2)=ones(size(i2)).*20e3*sqrt(2);
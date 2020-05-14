% Dispersion de vitesse des etoiles du disque dans la direction radiale
%
% Parametres : coordonnees Galactiques en pc
% Sortie: dispersion de vitesse en m/s

function res = siglp(R,z,t)
%res = ones(size(R)).*20e3*sqrt(2);
%res = ones(size(R)).*20e3;

global Rcoro Ro


% donnée de Pasetto, S. et al 
%"Thin disk kinematics from RAVE and the solar motion", Astronomy and Astrophysics 547 (2012), pp. A71.

res = ones(size(R)) * 27.4


\sigma_{r}^{thin}=27.4\pm1.1\ km/s\quad	\sigma_{r}^{thick}=56.1\pm3.8\ km/s
\sigma_{\theta}^{thin}=20.8\pm1.2\ km/s\quad	\sigma_{\theta}^{thick}=46.1\pm6.7\ km/s
\sigma_{z}^{thin}=16.3\pm2.2\ km/s\quad	\sigma_{z}^{thick}=35.1\pm3.4\ km/s


% i0 = find( R <= Rcoro );   
% i1 = find( R > Rcoro & R <= (Rcoro+Ro)/2 );
% i2 = find( R > (Rcoro+Ro)/2 );

% p=(20e3-110e3)*sqrt(2)/(Ro-Rcoro);
% res(i1)=p.*(R(i1)-Rcoro)+110e3*sqrt(2);

% res(i2)=ones(size(i2)).*20e3*sqrt(2);



%p=(20e3-110e3)*sqrt(2)/(Ro-Rcoro);
%l=length(R);
%for i=1:l;
%   if (R(i) <= Rcoro)
%    
%    elseif (R(i) > Rcoro & R(i) <= (Rcoro+Ro/2) );
%    res(i)= p.*(R(i)-Rcoro)+110e3*sqrt(2);
%    elseif (R(i) > (Rcoro+Ro)/2 );
%    res(i)=20e3*sqrt(2);
%    else disp(['erreur2']);
%    end;
%end;





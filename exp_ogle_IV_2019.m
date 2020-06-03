%---------------------------------------------------------------------------------------------------------------
%Microlensing optical depth and event rate in the OGLE-IV Galactic plane fields
%
% P. Mroz, A. Udalski, M.K. Szymanski, I. Soszynski, P. Pietrukowicz,
% S. Kozlowski, J. Skowron, R. Poleski, K. Ulaczyk, M. Gromadzki,
% K. Rybicki, P. Iwanek, and M. Wrona
%---------------------------------------------------------------------------------------------------------------

% Exposure ?


ttobs=tau/(gam/1e6/365.25);
disp(['<tobs> (en jours) = ' num2str(ttobs)]);
% N=gam*exposure;
% disp(['nb d''evt  = ' num2str(N)]);

taur=gam*pi/2*uT*mean(te)/365.25/1e6;
taur=real(taur);
disp(['tau (avec gamma integré par MC) = ' num2str(taur)]);


%-----------------
%Evenements
%-----------------

load ../OGLEIV/table3.dat

ogle_field = table_B1(:,2);
glon = table_B1(:,12)
glat = table_B1(:,13);
teff = table_B1(:,15);


%---------------
%Efficacité
%-----------
load ../OGLEIV/eff
%-----------------------------------------------------------------------------------------------
% Interpolation lineaire de l'efficacite pour determiner la probabilite qu'un evt a d'etre garde
%-----------------------------------------------------------------------------------------------

[tinterpmacho,indice] = sort(tmacho2005); 
effinterpmacho = effmacho2005t(indice) ;

%Il y a des doublons dans les données, il faut les supprimer pour que l'interpolation se passe correctement
indices = [1:length(tinterpmacho)-1];
il = find(tinterpmacho(indices)~=tinterpmacho(indices+1));
tinterpmacho = [tinterpmacho(il),tinterpmacho(length(tinterpmacho))];
effinterpmacho = [effinterpmacho(il),effinterpmacho(length(tinterpmacho))];


teffmaxm=max(tinterpmacho);
teffminm=min(tinterpmacho);

i1 = find((te<=teffmaxm)&(te>=teffminm));
effsimmacho = zeros(1,length(te));	% applique une efficacite nulle aux durees superieures et inferieures
effsimmachoblend = zeros(1,length(te));
effsimmacho(i1) = interp1(tinterpmacho,effinterpmacho,te(i1));
effsimmachoblend(i1) = interp1(tinterpmacho,effinterpmacho,teblend(i1));

%--------------------------------------------------------------------------------------------------------------------------
% compare le nombre aleatoire precedent a l'efficacite que l'on vient de calculer afin de decider si l'evt est garde ou non
%--------------------------------------------------------------------------------------------------------------------------

%tirage au sort pour l'efficacité
ramacho = rand(1,length(te))*max(effinterpmacho);

% On choisit l'efficacité ici en prenant les bons indices i

i = find(ramacho-effsimmacho<=0); 

teobs = te(i);

ib = find(ramacho-effsimmachoblend<=0); 

teobsblend = teblend(ib); % On récupère les éléments qui sont soumis au blending avec le calcul d'avant


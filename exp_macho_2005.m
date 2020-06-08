%---------------------------------------------------------------------------------------------------------------
% MACHO 8 Mai 2005 (accepté) : Microlensing optical depth toward the galactic bulge... Popowski & al
%---------------------------------------------------------------------------------------------------------------

% Exposure Macho 2005
exposure = 2530/365.25 * 1260000/10^6;


% ttobs=tau/(gam/1e6/365.25);
% disp(['<tobs> (en jours) = ' num2str(ttobs)]);
% N=gam*exposure;
% disp(['nb d''evt  = ' num2str(N)]);
% taur=gam*pi/2*uT*mean(te)/365.25/1e6;
% taur=real(taur);
% disp(['tau (avec gamma integré par MC) = ' num2str(taur)]);


%-----------------
%Evenements
%-----------------

%attention, les machos prennent le diametre d'einstein --> diviser les tmacho par 2 !

tmacho2005 = [165,153,107,38.3,287,45.6,14.2,254,53.5,30.8,152,17.5,45.6,47.3,45.8,66.9,9.26,18.5,28.7,28.7,19.8,27.9,6.53,36.9,23.7,64.9,14.7,25.5,9.77,53.7,23.8,126,13.6,74.4,19.4,41.7,60.9,30.8,24.3,13.7,312,44.7,60.6,72,125,20.1,74.6,5.67,13.7,5.46,50.1,19.3,47.6,115,12.3,54.3,15.9,25.5,48.5,11.6,9.59,15];
tmacho2005 = tmacho2005/2;
Amaxmacho2005 = [4.80,3.88,2.23,1.62,10.1,2.09,1.58,5.73,1.76,1.65,1.82,13.2,2.05,4.05,2.16,2.39,1.78,4.09,1.62,1.67,2.54,1.77,1.62,1.97,3.57,1.55,4.47,1.58,12.6,2.23,3.6,8.35,6.28,2.14,1.58,2.31,1.58,2.77,1.77,2.18,3.42,1.92,3.68,2.65,2.54,1.82,3.07,6.01,2.48,2.8,8.64,22.7,4.95,2.26,2.14,7.52,3.94,2.08,2.43,3.60,4.76,2.03];

effmacho2005tA = [0.68,0.39,0.55,0.48,0.77,0.46,0.41,0.78,0.48,0.47,0.66,0.40,0.56,0.57,0.57,0.58,0.41,0.53,0.44,0.45,0.41,0.53,0.29,0.59,0.44,0.46,0.45,0.54,0.45,0.56,0.54,0.73,0.49,0.56,0.27,0.22,0.33,0.30,0.26,0.25,0.58,0.29,0.31,0.34,0.38,0.24,0.32,0.098,0.26,0.13,0.32,0.29,0.31,0.38,0.24,0.34,0.29,0.3,0.31,0.21,0.2,0.25];
effmacho2005t = [0.55,0.32,0.43,0.38,0.63,0.39,0.33,0.62,0.39,0.4,0.52,0.31,0.48,0.47,0.48,0.46,0.34,0.42,0.4,0.4,0.33,0.46,0.27,0.49,0.36,0.39,0.32,0.44,0.34,0.46,0.43,0.58,0.39,0.48,0.21,0.2,0.24,0.24,0.24,0.19,0.44,0.24,0.25,0.25,0.33,0.20,0.25,0.059,0.2,0.086,0.26,0.24,0.26,0.32,0.20,0.26,0.21,0.24,0.26,0.16,0.14,0.20];

% Résultat du même article corrigé par les auteurs des effets du blending : attention forte barres d'erreur pour certains, j'ai enlevé les 3 événements super longs

temachoblending2005 = [196.6,148,103,46,68.7,14.2,43,78,133,17.84,46.2,31,63.3,12.7,16.7,3.92,25.7,16.5,19,30,24,122,16.8,15.7,93,26,121.3,12.48,88.5,15,48.7,61,30.6,3.68,13.1,454,38,61,83,48.5,7.5,13.4,11,64.7,21.8,80,20.7,119,45.3,10.9,14.5,13.9];
temachoblending2005 = temachoblending2005 / 2;

%Données
teff = tmacho2005;
eff = effmacho2005t;

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
ib = find((teblend<=teffmaxm)&(teblend>=teffminm));
effsimmacho = zeros(1,length(te));	% applique une efficacite nulle aux durees superieures et inferieures
effsimmachoblend = zeros(1,length(teblend));
effsimmacho(i1) = interp1(tinterpmacho,effinterpmacho,te(i1));
effsimmachoblend(ib) = interp1(tinterpmacho,effinterpmacho,teblend(ib));


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


%Script pour tracer des histogrammes des fonctions blendées à partir d'un modèle stocké dans un fichier txt classique 

%----------------------------------------
%recuperation des evenements selectionnes
%----------------------------------------

%On choisit le modèle de base

load ../graph/evenements.txt
x=evenements(:,1);
ds=evenements(:,2);
v=evenements(:,3);
m=evenements(:,4);
te=evenements(:,5);

x=x';
ds=ds';
v=v';
m=m';
te=te';


%---------------
%calcul de gamma
%---------------

gamma=tau/uT*2/pi/mean(te)*1e6*365.25;
disp(['gamma (calcule par le te moyen) = ' num2str(gamma)]);


gam1=4*sqrt(GMsol)/c*uT/sqrt(pc*pc*pc)*length(te)/(n*nbsimul)*86400*365.25*1e6;
gam=gam1*Gammax;
disp(['gamma (integre par MC) = ' num2str(gam)]);


ttobs=tau/(gam/1e6/365.25);
disp(['<tobs> (en jours) = ' num2str(ttobs)]);
% N=gam*exposure;
% disp(['nb d''evt  = ' num2str(N)]);

taur=gam*pi/2*uT*mean(te)/365.25/1e6;
taur=real(taur);
disp(['tau (avec gamma integré par MC) = ' num2str(taur)]);

temax = 30;
nbre_bin = temax;

%model
[hist_model, edges] = histcounts(teff, nbre_bin, 'BinLimits',[0,temax], 'Normalization', 'probability');

%------------------------
% Application du blending
%------------------------


% Param�tre pour le blending 
f = 0.5; % fraction des évenements concernés par le blending
nbar = 1.257; % P(n) = fonction(nbar) = f avec P(n) la proba d'avoir n étoiles dans DeltaS
script_blending 
[hist_05, edges] = histcounts(teff, nbre_bin, 'BinLimits',[0,temax], 'Normalization', 'probability');

f = 0.4
nbar = 1.6190
script_blending 
[hist_04, edges] = histcounts(teff, nbre_bin, 'BinLimits',[0,temax], 'Normalization', 'probability');

f = 0.6
nbar = 0.9474
script_blending 
[hist_06, edges] = histcounts(teff, nbre_bin, 'BinLimits',[0,temax], 'Normalization', 'probability');

for j =1:length(centre);
centre(j)=(edges(j)+edges(j+1))/2;
end

%plot 
figure(1)
hold on;
plot(centre, hist_model);
plot(centre, hist_04);
plot(centre, hist_05);
plot(centre, hist_06);

xlabel('t_{e}')
ylabel('Nombre d''évènements par unité de t_{e}')

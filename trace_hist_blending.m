%Script pour tracer des histogrammes des fonctions blendées à partir d'un modèle stocké dans un fichier txt classique 
clear;close all;
%----------------------------------------
%recuperation des evenements selectionnes
%----------------------------------------

%On choisit le modèle de base

load ../graph/evenements_1.txt

te=evenements_1(:,5);

te=te';


% %---------------
% %calcul de gamma
% %---------------
% 
% gamma=tau/uT*2/pi/mean(te)*1e6*365.25;
% disp(['gamma (calcule par le te moyen) = ' num2str(gamma)]);
% 
% 
% gam1=4*sqrt(GMsol)/c*uT/sqrt(pc*pc*pc)*length(te)/(n*nbsimul)*86400*365.25*1e6;
% gam=gam1*Gammax;
% disp(['gamma (integre par MC) = ' num2str(gam)]);
% 
% 
% ttobs=tau/(gam/1e6/365.25);
% disp(['<tobs> (en jours) = ' num2str(ttobs)]);
% % N=gam*exposure;
% % disp(['nb d''evt  = ' num2str(N)]);
% 
% taur=gam*pi/2*uT*mean(te)/365.25/1e6;
% taur=real(taur);
% disp(['tau (avec gamma integré par MC) = ' num2str(taur)]);

temax = 30;
nbre_bin = temax;

%model
[hist_model, edges] = histcounts(te, nbre_bin, 'BinLimits',[0,temax], 'Normalization', 'probability');
disp(['mean te = ' num2str(mean(te))])

%------------------------
% Application du blending
%------------------------

global Vinf Vsup normfl uT

uT = 1;
Vinf = 22; % magnitudes limites des étoiles observées.
Vsup = 16;

Linf=lum(Vinf);
Lsup=lum(Vsup);

% Norme de la fonction de luminosité :

normfl = integral(@fl, Linf, Lsup);

% Param�tre pour le blending 
f = 0.5; % fraction des évenements concernés par le blending
nbar = 1.257; % P(n) = fonction(nbar) = f avec P(n) la proba d'avoir n étoiles dans DeltaS
script_blending 
[hist_05, edges] = histcounts(teblend, nbre_bin, 'BinLimits',[0,temax], 'Normalization', 'probability');

disp(['mean te = ' num2str(mean(teblend))])

f = 0.4;
nbar = 1.6190;
script_blending 
[hist_04, edges] = histcounts(teblend, nbre_bin, 'BinLimits',[0,temax], 'Normalization', 'probability');
disp(['mean te = ' num2str(mean(teblend))])

f = 0.6;
nbar = 0.9474;
script_blending 
[hist_06, edges] = histcounts(teblend, nbre_bin, 'BinLimits',[0,temax], 'Normalization', 'probability');
disp(['mean te = ' num2str(mean(teblend))])

centre = zeros(size(edges)-[0,1]);
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
legend('modèle', 'f = 0.4', 'f = 0.5', 'f = 0.6')
xlabel('t_{e}')
ylabel('Nombre d''évènements par unité de t_{e}')

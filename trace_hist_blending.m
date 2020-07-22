%Script pour tracer des histogrammes des fonctions blendées à partir d'un modèle stocké dans un fichier txt classique 
clear;close all;
%----------------------------------
% Constantes physiques (unites SI)
%----------------------------------
G=6.672e-11;	pc=3.08567802e16;
kpc=pc*1e3;  	Msol=1.989e30;
c=299792458;	GMsol=1.32712497e20;

%----------------------------------------
%recuperation des evenements selectionnes
%----------------------------------------

%On choisit le modèle de base

load ../graph/modif_FM/evenements_1.txt

te=evenements_1(:,5);

te=te';

fid = fopen('../graph/modif_FM/simul_para_1.txt');
tau = str2double(fgets(fid));
n = str2double(fgets(fid));
nbsimul = str2double(fgets(fid));
tau_1 = str2double(fgets(fid));
Gammax = str2double(fgets(fid));
uT = str2double(fgets(fid));
At = str2double(fgets(fid));
mmean = str2double(fgets(fid));

%---------------
%calcul de gamma
%---------------

gamma=tau/uT*2/pi/mean(te)*1e6*365.25;
disp(['gamma (calcule par le te moyen) = ' num2str(gamma)]);


gam1=4*sqrt(GMsol)/c*uT/sqrt(pc*pc*pc)*length(te)/(n*nbsimul)*86400*365.25*1e6/mmean;
gam=gam1*Gammax;
disp(['gamma (integre par MC) = ' num2str(gam)]);


ttobs=tau/(gam/1e6/365.25);
disp(['<tobs> (en jours) = ' num2str(ttobs)]);
% N=gam*exposure;
% disp(['nb d''evt  = ' num2str(N)]);

taur=gam*pi/2*uT*mean(te)/365.25/1e6;
taur=real(taur);
disp(['tau (avec gamma integré par MC) = ' num2str(taur)]);

%histogrammes
temax = 100;
nbre_bin = temax/2;

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
[hist_05, edges] = histcounts(teblend(teblend ~=0), nbre_bin, 'BinLimits',[0,temax], 'Normalization', 'probability');

disp(['mean te = ' num2str(mean(teblend))])

f = 0.4;
nbar = 1.6190;
script_blending 
[hist_04, edges] = histcounts(teblend(teblend ~=0), nbre_bin, 'BinLimits',[0,temax], 'Normalization', 'probability');
disp(['mean te = ' num2str(mean(teblend))])

f = 1;
nbar = 0;
script_blending 
[hist_06, edges] = histcounts(teblend(teblend ~=0), nbre_bin, 'BinLimits',[0,temax], 'Normalization', 'probability');
disp(['mean te = ' num2str(mean(teblend))])
gam1=4*sqrt(GMsol)/c*uT/sqrt(pc*pc*pc)*length(teblend)/(n*nbsimul)*86400*365.25*1e6/mmean;
gam1*Gammax
gam*pi/2*uT*mean(teblend)/365.25/1e6

%%
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

%---------------
%graph pour tau
%-----------------

%fraction en fonction de n :
fraction = @(X) X.*exp(-X)./(1-exp(-X));

n_frac = (0:1e-3:1)*5;
tau_b = zeros(size(n_frac));
gam_b = zeros(size(n_frac));
fract = fraction(n_frac);

for i = 1:length(n_frac)
    f = fract(i); nbar = n_frac(i);
    script_blending
    
    gam1=4*sqrt(GMsol)/c*uT/sqrt(pc*pc*pc)*length(teblend)/(n*nbsimul)*86400*365.25*1e6/mmean;
    gam_b(i)=gam1*Gammax;

    tau_b(i)=gam*pi/2*uT*mean(teblend)/365.25/1e6;
end

figure(2)
plot(fract, tau_b/taur)
% Programme de micro-lentilles gravitationnelles (bulbe et disque)

clear


global vlimit

vlimit = 1000e3;

%----------------------
% Nombre de simulations
%----------------------

n = 60e5;
% n = 5000;
nbsimul=10; %a augmenter pour meilleure stat

%----------------------------------------------------------------------
% Param�tres de la fonction de distribution de la distance de la source
%----------------------------------------------------------------------

global dsup dinf 

dsup = 20000.;
dinf = 800.;
%distance en parsec

%-------
% soleil 
%-------
global Ro elev
Ro = 8000;	   % distance Sun-GC en pc
elev = 26 ;      % elevation au dessus du plan du disque
Lkpc=Ro/1000;  % distance Sun-GC en kpc


%----------------------------
%rayon de corotation (en pc)
%----------------------------

global Rcoro
Rcoro = 3500;


%---------------------------------------       
% param�tres du calcul de microlentilles
%---------------------------------------
global l b

% definition de la fenetre de Baade dans la majeure partie des articles :  l = 1 et b = -4
% A priori c'est cette definition qui est juste.
l = 1 *pi/180;    % direction d'observation en radian
b = -4 *pi/180;


% definition de la fenetre de Baade dans les theses de Mera et Alibert : l = 4 et b = -1
% l = 4 *pi/180;    % direction d'observation en radian
% b = -1 *pi/180;

uT = 1;		   % Seuil de d�tection en param�tre d'impact
AT = 3/sqrt(5);    % Seuil de d�tectiocreen en amplification


%-----------------------------------
% Param�tres de la fonction de masse
%-----------------------------------
global minf msup

minf=0.01;
msup=100;

%----------------------------------------
% Param�tres de la fonction de luminosité pour le blending (Alibert et al.)
%----------------------------------------

global Vinf Vsup 

Vinf = 22; % magnitudes limites des étoiles observées.
Vsup = 16;



%-------------------
% Monte-Carlo
%-------------------
tic
main
toc
%%
%----------------------------------------
% recuperation des evenements selectionnes
%----------------------------------------

load evenements.txt
x=evenements(:,1);
ds=evenements(:,2);
v=evenements(:,3);
m=evenements(:,4);
te=evenements(:,5);

te=te';
x = x';
ds = ds';

load ../graph/modif_FM/evenements_1.txt
% x=evenements_1(:,1);
% ds=evenements_1(:,2);
% v=evenements_1(:,3);
% m=evenements_1(:,4);
te_model=evenements_1(:,5);


te_model = te_model';
fid = fopen('../graph/modif_FM/simul_para_1.txt');
tau = str2double(fgets(fid));
n = str2double(fgets(fid));
nbsimul = str2double(fgets(fid));
tau_1 = str2double(fgets(fid));
Gammax = str2double(fgets(fid));
uT = str2double(fgets(fid));
At = str2double(fgets(fid));
mmean = str2double(fgets(fid));


disp(' ')
close all
%---------------
%calcul de gamma
%---------------

gamma=tau/uT*2/pi/mean(te)*1e6*365.25;
disp(['gamma (calcule par le te moyen) =    ' num2str(gamma)]);


gam1=4*sqrt(GMsol)/c*uT/sqrt(pc*pc*pc)*length(te)/(n*nbsimul)*86400*365.25*1e6;
gam=gam1*Gammax;
disp(['gamma (integre par MC) = ' num2str(gam)]);


ttobs=tau/(gam/1e6/365.25);
disp(['<tobs> (en jours) = ' num2str(ttobs)]);

taur=gam*pi/2*uT*mean(te)/365.25/1e6;
taur=real(taur);
disp(['tau (avec gamma integré par MC) = ' num2str(taur)]);

%------------------------
% Application du blending
%------------------------

% Param�tre pour le blending 
f = 0.05;   %fraction des évenements unblendé f = P(1)
nbar = 4.51; % P(n) = fonction(nbar) = f avec P(n) la proba d'avoir n étoiles dans DeltaS

f = 0.5;
nbar = 1.257;
% f = 0.2;
% nbar = 2.6;
% f = 1;
% nbra = 0;
%retourn teblend (histogramme corrigé) et taurblend (profondeur optique corrigée)
script_blending 

%-----------
%Choix de l'expérience à analyser
%donne teff : tet des observations et eff : efficacité
%------------

exp_ogle_III_2015

%---------------
%calcul de gamma
%---------------

disp('Grandeurs avec intervention de l''efficacite experimentale :')

%------------------------------------------------------------------------------------------------
% on ne peut pas calculer le gamma par la formule avec le tobs, car le tau ne prend pas en compte
% l'efficacite. Par contre, on peut deduire tau experimental a partir du gamma calcule par MC
%------------------------------------------------------------------------------------------------

gamobs = gam/length(te)*length(teobs);
disp(['gamma (integre par MC) = ' num2str(gamobs)]);
% 
tauobs=gamobs*pi/2*uT*mean(teobs)/365.25/1e6;
disp(['tau obs (calcule par le te moyen) = ' num2str(tauobs)]);
% 
% tauobsblend=tauobs * gmean * (nbar/(1-exp(-nbar)));
% tauobsblend=real(tauobsblend);
% % disp(['tau observé avec blending (Alibert 2005)  = ' num2str(tauobsblend)]);
% 
gamobsb = gam/length(te)*length(teobsblend);
disp(['gamma avec blending (integre par MC) = ' num2str(gamobsb)]);
% 
tauobsb=gamobsb*pi/2*uT*mean(teobsblend)/365.25/1e6;
disp(['tau obs avec blending (calcule par le te moyen) = ' num2str(tauobsb)]);
% 
% disp(['rapport tau_blend/tau_obs_théorique = ' num2str(tauobsb/tauobs)]);

%------------------------
% affichage des resultats
%------------------------


%Paramètre graph
bin_max = 100;
nbre_bin = bin_max/2;


%Trace la distribde te  pour le modèle et la courbe stockée localement
[hist_1, edges] = histcounts(te, nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');
[hist_model, edges] = histcounts(te_model, nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');

%tracé distribution avec blending uniquement la courbe stockée localement
[histb, edges] = histcounts(teblend, nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');

%Courbe expérimentale (avec l'efficacité) :
[hist_obs, edges] = histcounts(teobs, nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');
[hist_obs_b, edges] = histcounts(teobsblend, nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');

%courbe de l'expérience
[hist_exp_err, edges] = histcounts(teff, nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');
[hist_exp_BW, edges] = histcounts(teff(i_BW), nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');
[hist_exp_fs, edges] = histcounts(teff(i_fs), nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');

centre = zeros(size(edges)-[0,1]);

for j =1:length(centre);
centre(j)=(edges(j)+edges(j+1))/2;
end

%Graph normalisé
figure(16)
hold on;
plot(centre, hist_1, 'black');
plot(centre, hist_model, 'red');
legend("local", "model")
title('comparaison local et modèle')
xlabel('t_{e}')
ylabel('Nombre d''évènements par unité de t_{e}')



%Graph normalisé expérience et exp simulée avec l'efficacité pour OGLE III
figure(18)
hold on;
plot(centre, hist_obs, 'red');
plot(centre, hist_obs_b, 'black');
M = length(hist_exp_err);
plot(edges(sort([1:M 1:M])), [0 , 0, hist_exp_err(sort([1:M 2:M-1]))])
M = length(hist_exp_BW);
plot(edges(sort([1:M 1:M])), [0 , 0, hist_exp_BW(sort([1:M 2:M-1]))], 'g')
M = length(hist_exp_fs);
plot(edges(sort([1:M 1:M])), [0 , 0, hist_exp_fs(sort([1:M 2:M-1]))])
legend('hist modèle', 'hist modèle avec blending (f=0.5)', 'OGLE III (all stars)', 'OGLE III (fenêtre de Baade)', 'OGLE III(non blendé, f_{s}>0.2)')
xlabel('t_{e}')
ylabel('Nombre d''évènements par unité de t_{e}')

%%
%Graph log
te_min = min(teobs); te_max = max(teobs); M = 30;
edges_log=te_min*(te_max/te_min).^([0:M]/M);
x=edges_log(sort([1:M+1 1:M+1])); 

hist_exp_log = histcounts(teff,edges_log);
hist_exp_log_fs = histcounts(teff(i_fs),edges_log);
hist_exp_log_BW = histcounts(teff(i_BW),edges_log);
hist_obs_log = histcounts(teobs,edges_log);
hist_obs_b_log = histcounts(teobsblend,edges_log);

centre = zeros(size(edges_log)-[0,1]);
for j =1:length(centre);
centre(j)=(edges_log(j)+edges_log(j+1))/2;
end

if ishandle(2)
    close(2)
end
figure(2)
hold on
plot(centre, hist_obs_log./length(teobs))
plot(centre, hist_obs_b_log./length(teobsblend))
plot(x, [0 hist_exp_log(sort([1:M 1:M])) 0]./length(teff))
% plot(x, [0; hist_exp_log_BW(sort([1:M-1 1:M-1])); 0]/length(teff(i_BW)))
plot(x, [0 hist_exp_log_fs(sort([1:M 1:M])) 0]/length(teff(i_fs)), 'b')
errorbar(centre, hist_exp_log_fs./length(teff(i_fs)), 1./sqrt(hist_exp_log_fs)./length(teff(i_fs)), 'b.')

legend('hist modèle', 'hist modèle avec blending (f=0.5)', 'OGLE III (all stars)', 'OGLE III(non blendé, f_{s}>0.2)')
legend('Location', 'best')
xlabel('log(t_{e})')
ylabel('Nombre d''évènements par unité de t_{e} (échelle log)')
% title('normalisé pour expo = 10 et gamobs')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
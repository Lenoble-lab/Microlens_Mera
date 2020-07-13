% Programme de micro-lentilles gravitationnelles (bulbe et disque)

clear


global vlimit

vlimit = 1000e3;

%----------------------
% Nombre de simulations
%----------------------

n = 20000;
% n = 5000;
nbsimul=500; %a augmenter pour meilleure stat

%----------------------------------------------------------------------
% Param�tres de la fonction de distribution de la distance de la source
%----------------------------------------------------------------------

global dsup dinf 

dsup = 15000.;
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
AT = 3/sqrt(5);    % Seuil de d�tection en amplification


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
main


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
star_pop = evenements(:,6);

x=x';
ds=ds';
v=v';
m=m';
te=te';
star_pop = star_pop';

disp(' ')
close all

%------------------
%Test pour la PDMF
%--------------

figure(40)
hold on
[Y, edges] = histcounts(m_tot_pdmf, 1000, 'BinLimits',[0,5]);
[Y_bon, edges] = histcounts(m, 1000, 'BinLimits',[0,5]);

M = length(Y);
plot(edges(sort([1:M 2:M+1])), Y(sort([1:M 1:M]))/length(m_tot_pdmf)/2.7./sqrt(edges(sort([1:M 1:M]))))
plot(edges(sort([1:M 2:M+1])), Y(sort([1:M 1:M]))/length(m_tot_pdmf))
plot(edges(sort([1:M 2:M+1])), Y_bon(sort([1:M 1:M]))/length(m))
legend('PDMF avec facteu', 'PDMF brut','eve retenu brut')
title("PDMF")
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')


% axis([0 10 1e-2 1e3])

disp(['fraction de rémanents en nombre (BD, MS, WD, NS, BH) ', num2str(mean(frac_N_tot))])
disp(['fraction de rémanents en masse (BD, MS, WD, NS, BH) ', num2str(mean(frac_M_tot))])
disp(['fraction d evenements (BD, MS, WD, NS, BH) ', num2str(mean(frac_eve_tot))])

%-------------
%test distrib vitesse
%------------------
% figure(41)
% hold on
% [v_bu_bu, edges] = histcounts(v_tot(intersect(ibul, ibus)), 500);
% [v_dm_bu, edges] = histcounts(v_tot(intersect([idml, idel], ibus)), 500);
% [v_dm_dm, edges] = histcounts(v_tot(intersect([idml, idel], [idms, ides])), 500);
% 
% v_dm = v_tot(intersect([idml, idel], [idms, ides]));
% v_test = v_dm(find(ds(intersect([idml, idel], [idms, ides]))>8e3));
% [v_dm_dm_1, edges_1] = histcounts(v_test, 500);
% 
% [v_retenu, edges] = histcounts(v, 500);
% 
% M = length(Y);
% plot(edges(sort([1:M 2:M+1])), v_bu_bu(sort([1:M 1:M])))
% plot(edges(sort([1:M 2:M+1])), v_dm_bu(sort([1:M 1:M])))
% plot(edges(sort([1:M 2:M+1])), v_dm_dm(sort([1:M 1:M])))
% % plot(edges_1(sort([1:M 2:M+1])), v_dm_dm_1(sort([1:M 1:M])))
% % plot(edges(sort([1:M 2:M+1])), v_retenu(sort([1:M 1:M])))
% legend('bulbe/bulbe', 'bulbe/disque', 'disque/disque')
% title("vitesse")
% set(gca, 'YScale', 'log')
% % axis([2e3 10e3 0 3e6])



%-------------
%test distrib L
%------------------
% figure(42)
% hold on
% [Y, edges] = histcounts(L_tot, 500);
% M = length(Y);
% plot(edges(sort([1:M 2:M+1])), Y(sort([1:M 1:M]))/length(m))
% title("Distance source")
% set(gca, 'YScale', 'log')

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


temax = 100;
nbre_bin = temax/5;

% exp_ogle_2006
% teff_ogle = teff;
% [hist_ogle, edges] = histcounts(teff, nbre_bin, 'BinLimits',[0,temax]);
% 
% 
% exp_macho_2005
% [hist_macho, edges] = histcounts(teff, nbre_bin, 'BinLimits',[0,temax]);
% sort(teff)
% figure(1)
% bar(edges(1:end-1),hist_macho)

% exp_eros_2006
% % teff_eros = teff;
% [hist_eros, edges] = histcounts(teff, nbre_bin, 'BinLimits',[0,temax]);
% 
% figure(1)
% bar(edges(1:end-1),[hist_ogle; hist_macho; hist_eros]')
% legend('OGLE', 'MACHO', 'EROS')

% exp_ogle_IV_2019
% exp_MOA_2016
exp_ogle_III_2015
% exp_ogle_II_2006

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


%telechargement de la courbe du modèle

load ../graph/07_07_evenements/evenements_1.txt
te_model = evenements_1(:,5);

%Paramètre graph
bin_max = 100;
nbre_bin = bin_max;


%Trace la distribde te  pour le modèle et la courbe stockée localement
[hist_local, edges] = histcounts(te, nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');
[hist_model, edges] = histcounts(te_model, nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');

%Selon la pop d'étoiles
[hist_BD, edges] = histcounts(te(star_pop == 1), nbre_bin, 'BinLimits',[0,bin_max]);
[hist_MS, edges] = histcounts(te(star_pop == 2), nbre_bin, 'BinLimits',[0,bin_max]);
[hist_WD, edges] = histcounts(te(star_pop == 3), nbre_bin, 'BinLimits',[0,bin_max]);
[hist_NS, edges] = histcounts(te(star_pop == 4), nbre_bin, 'BinLimits',[0,bin_max]);
[hist_BH, edges] = histcounts(te(star_pop == 5), nbre_bin, 'BinLimits',[0,bin_max]);



%tracé distribution avec blending uniquement la courbe stockée localement
[histb, edges] = histcounts(teblend, nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');

%Courbe expérimentale (avec l'efficacité) :
[hist_obs, edges] = histcounts(teobs, nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');
[hist_obs_b, edges] = histcounts(teobsblend, nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');

%courbe de l'expérience
[hist_exp_err, edges] = histcounts(teff, nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');
% [hist_exp_BW, edges] = histcounts(teff, nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');


centre = zeros(size(edges)-[0,1]);

for j =1:length(centre);
centre(j)=(edges(j)+edges(j+1))/2;
end

%Graph normalisé
figure(16)
hold on;
plot(centre/66, hist_local/max(hist_local), 'black');
plot(centre/66, hist_model/max(hist_model), 'red');
% plot(centre, hist_BD/length(find(star_pop == 1)));
% plot(centre, hist_MS/length(find(star_pop == 2)));
% plot(centre, hist_WD/length(find(star_pop == 3)));
% plot(centre, hist_NS/length(find(star_pop == 4)));
% plot(centre, hist_BH/length(find(star_pop == 5)));

title('comparaison local et modèle')
legend('local', 'model', 'BD', 'MS', 'WD', 'NS', 'BH')
xlabel('t_{e}')
ylabel('Nombre d''évènements par unité de t_{e}')
% set(gca, 'XScale', 'log')

%%
%graph en fonction de l'exposition
exposure = 10;
figure(17)
hold on;
plot(centre, hist_local/length(te)*gamobs*exposure, 'red');
% plot(centre, hist_obs_b*gamobsb*exposure, 'black');
% plot(centre, (hist_BD+hist_WD+hist_MS)/length(te)*gamobs*exposure);

plot(centre, hist_BD/length(te)*gamobs*exposure);
plot(centre, hist_MS/length(te)*gamobs*exposure);
plot(centre, hist_WD/length(te)*gamobs*exposure);
plot(centre, hist_NS/length(te)*gamobs*exposure);
plot(centre, hist_BH/length(te)*gamobs*exposure);
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
legend('hist modèle', 'BD', 'MS', 'WD', 'NS', 'BH')
xlabel('t_{e}')
ylabel('Nombre d''évènements par unité de t_{e}')

%graph logx
figure(22)
hold on
M = 50;
vmin = 0.1; vmax = 700;
edges_log=vmin*(vmax/vmin).^([0:M]/M);
x=edges_log(sort([1:M 1:M])); 

log_BD = histc(te(star_pop == 1),edges_log);
log_MS = histc(te(star_pop == 2),edges_log);
log_WD = histc(te(star_pop == 3),edges_log);
log_NS = histc(te(star_pop == 4),edges_log);
log_BH = histc(te(star_pop == 5),edges_log);
log_tot = histc(te,edges_log);

% log_BD = histc(intersect(te(star_pop == 1), i_eff),edges_log);
% log_MS = histc(intersect(te(star_pop == 2), i_eff),edges_log);
% log_WD = histc(intersect(te(star_pop == 3), i_eff),edges_log);
% log_NS = histc(intersect(te(star_pop == 4), i_eff),edges_log);
% log_BH = histc(intersect(te(star_pop == 5), i_eff),edges_log);
% log_tot = histc(te(i_eff),edges_log);

plot(x, [0 log_BD(sort([1:M-1 1:M-1])) 0]/length(te)*gamobs*exposure)
plot(x, [0 log_MS(sort([1:M-1 1:M-1])) 0]/length(te)*gamobs*exposure)
plot(x, [0 log_WD(sort([1:M-1 1:M-1])) 0]/length(te)*gamobs*exposure)
plot(x, [0 log_NS(sort([1:M-1 1:M-1])) 0]/length(te)*gamobs*exposure)
plot(x, [0 log_BH(sort([1:M-1 1:M-1])) 0]/length(te)*gamobs*exposure)
plot(x, [0 log_tot(sort([1:M-1 1:M-1])) 0]/length(te)*gamobs*exposure)


legend('BD', 'MS', 'WD', 'NS', 'BH', 'tot', 'expérience')
xlabel('log(t_{e})')
ylabel('Nombre d''évènements par unité de t_{e} (échelle log)')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')


%expérience
figure(18)
hold on;
% plot(centre, hist_obs, 'red');
% plot(centre, hist_obs_b, 'black');
log_BD_eff = histc(te(intersect(find(te(star_pop == 1)), i_eff)),edges_log);
log_MS_eff = histc(te(intersect(find(te(star_pop == 2)), i_eff)),edges_log);
log_WD_eff = histc(te(intersect(find(te(star_pop == 3)), i_eff)),edges_log);
log_NS_eff = histc(te(intersect(find(te(star_pop == 4)), i_eff)),edges_log);
log_BH_eff = histc(te(intersect(find(te(star_pop == 5)), i_eff)),edges_log);
log_tot_eff = histc(te(i_eff),edges_log);
log_exp_eff = histc(teff, edges_log);

plot(x, [0 log_BD_eff(sort([1:M-1 1:M-1])) 0]/length(te(i_eff)))
plot(x, [0 log_MS_eff(sort([1:M-1 1:M-1])) 0]/length(te(i_eff)))
plot(x, [0 log_WD_eff(sort([1:M-1 1:M-1])) 0]/length(te(i_eff)))
plot(x, [0 log_NS_eff(sort([1:M-1 1:M-1])) 0]/length(te(i_eff)))
plot(x, [0 log_BH_eff(sort([1:M-1 1:M-1])) 0]/length(te(i_eff)))
plot(x, [0 log_tot_eff(sort([1:M-1 1:M-1])) 0]/length(te(i_eff)))
plot(x, [0; log_exp_eff(sort([1:M-1 1:M-1])); 0]/length(teff))


legend('BD', 'MS', 'WD', 'NS', 'BH', 'tot', 'expérience')
xlabel('log(t_{e})')
ylabel('Nombre d''évènements par unité de t_{e} (échelle log)')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

% legend('hist modèle', 'hist modèle avec blending (f=0.5)', 'expérience')
% xlabel('t_{e}')
% ylabel('Nombre d''évènements par unité de t_{e}')



% temaxgraphe=100;  % pour le trace des graphe
% temax=5000;       % pour le calcul integral
% 
% % Definition des differents nombres de bin pour les graphiques
% nbbin1=temaxgraphe;
% nbbin2=temaxgraphe/2;     
% nbbin5=temaxgraphe/5;
% 
% bords1 = zeros(1, nbbin1) ; centre1 = zeros(1, nbbin1-1);
% bords2 = zeros(1, nbbin2) ; centre2 = zeros(1, nbbin2-1);
% bords5 = zeros(1, nbbin5) ; centre5 = zeros(1, nbbin5-1);
% 
% teobs = te;
% %Definition des tableaux bords et centre pour une largeur de 1 jour
% for j =1:nbbin1;
% bords1(j)=temaxgraphe/nbbin1*(j-1);
% end;
% for j =1:nbbin1-1;
% centre1(j)=(bords1(j)+bords1(j+1))/2;
% end
% centre1(nbbin1)=bords1(nbbin1)+(bords1(nbbin1)-bords1(nbbin1-1))/2;
% 
% i=find(teobs<temaxgraphe); % ce tableau i est le meme pour les trois graphes
% 
% h1=histc(teobs(i),bords1)*gamobs*exposure/length(teobs(i));
% 
% %Definition des tableaux bords et centre pour une largeur de 2 jour
% for j =1:nbbin2;
% bords2(j)=temaxgraphe/nbbin2*(j-1);
% end;
% for j =1:nbbin2-1;
% centre2(j)=(bords2(j)+bords2(j+1))/2;
% end
% centre2(nbbin2)=bords2(nbbin2)+(bords2(nbbin2)-bords2(nbbin2-1))/2;
% 
% h2=histc(teobs(i),bords2)*gamobs*exposure/length(teobs(i));
% 
% %Definition des tableaux bords et centre pour une largeur de 5 jour
% for j =1:nbbin5;
% bords5(j)=temaxgraphe/nbbin5*(j-1);
% end;
% for j =1:nbbin5-1;
% centre5(j)=(bords5(j)+bords5(j+1))/2;
% end
% centre5(nbbin5)=bords5(nbbin5)+(bords5(nbbin5)-bords5(nbbin5-1))/2;
% 
% h5=histc(teobs(i),bords5)*gamobs*exposure/length(teobs(i));
% 
% 
% % Tracage de ces trois graphes ensemble
% figure(18);
% hold off;
% hold on;
% plot(centre1,h1,'b-');
% plot(centre2,h2,'r-');
% plot(centre5,h5,'g-');
% axis([0,30,0,20]);
% title('Diagramme des evenements obtenu par simulation pour des largeurs de representation differentes');
% legend('largeur de 1 jour','largeur de 2 jours','largeur de 5 jours');
% xlabel('durees d''evenements');
% ylabel('nombre d''evenements par unite de temps');
% hold off;
% 
% 
% % TRACAGE des trois graphiques precedents separement avec la distribution experimentale correspondante
% 
% figure(19);
% j=find(teff<temaxgraphe);  % c'est le meme tableau pour les trois graphes
% hist(teff(j),bords1);% essayer aussi hi1=hist(teff(i),bords1)*gamobs*exposure/length(teff(i));
% teff(j);
% %bar(centre1,hi1,'r');
% hold on;
% plot(centre1,h1,'r');
% title('comparaison experience/simulation pour une largeur de 1 jours');
% legend('experience','simulation');
% hold off;
% 
% figure(20);
% hist(teff(j),bords2);   % essayer aussi hi2=hist(teff(i),bords2)*gamobs*exposure/length(teff(i));
% %bar(centre2,hi2,'r');
% hold on;
% plot(centre2,h2,'r');
% title('comparaison experience/simulation pour une largeur de 2 jours');
% legend('experience','simulation');
% hold off;
% 
% figure(21);
% hist(teff(j),bords5);   % essayer aussi hi5=hist(teff(i),bords5)*gamobs*exposure/length(teff(i));
% %bar(centre5,hi5,'r');
% hold on;
% plot(centre5,h5,'r');
% title('comparaison experience/simulation pour une largeur de 5 jours');
% legend('experience','simulation');
% hold off;
% 
% 
% 
% 
% % Calcul du nombre total d'evenement pour les differents cas
% 
% %calcul du nombre total d'evenements a partir de gamobs et l'exposure
% nombre1=exposure*gamobs;
% %calcul du nombre total d'evenements a partir de la simulation pour une largeur de 1 jour
% nbbin1=temax;
% for j =1:nbbin1;
% bords1(j)=temax/nbbin1*(j-1);
% end;
% i=find(teobs<temax); %c'est le meme tableau pour les trois (quatres) largeurs
% hn1=histc(teobs(i),bords1)*gamobs*exposure/length(teobs(i));
% diff=length(teobs)-length(i); % c'est la meme valeur pour tous
% nombre2=sum(hn1);
% %calcul du nombre total d'evenements a partir de la simulation pour une largeur de 2 jours
% nbbin2=temax/2;     
% for j =1:nbbin2;
% bords2(j)=temax/nbbin2*(j-1);
% end;
% hn2=histc(teobs(i),bords2)*gamobs*exposure/length(teobs(i));
% nombre3=sum(hn2);
% %calcul du nombre total d'evenements a partir de la simulation pour une largeur de 5 jours     
% nbbin5=temax/5;
% for j =1:nbbin5;
% bords5(j)=temax/nbbin5*(j-1);
% end;
% hn5=histc(teobs(i),bords5)*gamobs*exposure/length(teobs(i));
% nombre4=sum(hn5);
% %calcul du nombre total d'evenements a partir de la simulation pour une largeur de 0.5 jour
% nbbin=2*temax;
% for j =1:nbbin;
% bords(j)=temax/nbbin*(j-1);
% end;
% hn=histc(teobs(i),bords)*gamobs*exposure/length(teobs(i));
% nombre5=sum(hn);
% %calcul du nombre total d'evenements a partir des donnes experimentales
% nombre6=length(teff);
% 
% %Affichage de ces resultats
% disp(['Il n''y a pas eu prise en compte de ' num2str(diff) 'elements']);
% disp(['nb d''evt avec gamobs :' num2str(nombre1) ]);
% disp(['nb d''evt avec 1 j :' num2str(nombre2) ]);
% disp(['nb d''evt avec 2 j :' num2str(nombre3) ]);
% disp(['nb d''evt avec 5 j :' num2str(nombre4) ]);
% disp(['nb d''evt avec 0.5 j :' num2str(nombre5) ]);
% disp(['nb d''evt par l''experience :' num2str(nombre6) ]);
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Trace dgamma (c'est le dgamma global et non pas celui apres le Monte-Carlo)
% % Pour avoir quelque chose de coherent, trace plutot dgammaccepte
% % figure(22);
% % mgam=max(dgamma);
% % y=0:mgam/100000:mgam;
% % h4=histc(dgamma,y);
% % h4=h4+1;
% % h4=log(h4);
% % bar(y,h4);
% % title('distrib de dgamma');
% % hold off;
% %trace dgamma avec une echelle log pour la lisibilite
% %figure(23);
% %%dgamma=dgamma*1000+1;
% %ldg=log(dgamma)/log(1000);
% %hist(ldg,100000);
% %title('distrib de dgamma echelle log1000');
% %hold off;
% %n'a pas trop d'interet
% 
% 
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Densite (Avec les espacements pour ramener les dimensions)
% 
% %figure(24);
% %l=length(teff);
% %milieu(1)=0;
% %for i=2:l
% %%    milieu(i)=(teff(i-1)+teff(i))/2;
% %end;
% %for i=1:(l-1)
% %    ecart(i)=milieu(i+1)-milieu(i);
% %end;
% %ecart(l)=(teff(l)-milieu(l))*2;
% %h3=1./ecart;
% %bar(t,h3,'g');
% %title('idee bizarre');
% %hold off;
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % En essayant de comparer des "resultats" normalises
% % Ne sert a rien s'il n'y a pas une coherence avant
% %figure(25);
% %m1=max(h);
% %h1=h1/m1;
% %plot(centre,h1,'b-');
% %hold on;
% %h2=histc(teff,bords);
% %m2=max(h2);
% %h2=h2/m2;
% %plot(centre,h2,'r-.');  ne pas tracer avec celui la
% %bar(centre,h2,'g');
% %title('comparaison en normalisee');
% %hold off;
% 

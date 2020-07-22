% Programme de micro-lentilles gravitationnelles (bulbe et disque)

clearvars -except field field_li

%---------------------------------------------------------------------------------------------------------------
%Microlensing optical depth and event rate in the OGLE-IV Galactic plane fields
%
% P. Mroz, A. Udalski, M.K. Szymanski, I. Soszynski, P. Pietrukowicz,
% S. Kozlowski, J. Skowron, R. Poleski, K. Ulaczyk, M. Gromadzki,
% K. Rybicki, P. Iwanek, and M. Wrona
%---------------------------------------------------------------------------------------------------------------

%-----------------------------------------------
%Table 6. Basic information about analyzed fields
%------------------------------------------------
delimiter = ' ';
VarNames_table6 = {'field', 'ra', 'dec', 'glon', 'glat', 'N_stars', 'N_epochs'};
VarTypes_table6 = {'string', 'double', 'double', 'double', 'double', 'double', 'double'}; 

opts = delimitedTextImportOptions('VariableNames',VarNames_table6,'VariableTypes',VarTypes_table6,...
                                'Delimiter',delimiter, 'DataLines', 18, ...
                       'WhiteSpace', ' ', 'ConsecutiveDelimitersRule', 'join');
table6 = readtable('../OGLEIV/table6.dat',opts);
%-----------------------
%Table 4 Number of events detected in individual timescale bins.
%------------------------------

delimiter = ' ';
VarNames_table4 = {'bin', 'log_tE', 'BLG500', 'BLG501', 'BLG504', 'BLG505', 'BLG506', 'BLG511', 'BLG512', 'BLG534', 'BLG611'};
VarTypes_table4 = {'double', 'double', 'double', 'double', 'double', 'double', 'double',  'double', 'double', 'double', 'double'}; 

opts = delimitedTextImportOptions('VariableNames',VarNames_table4,'VariableTypes',VarTypes_table4,...
                                'Delimiter',delimiter, 'DataLines', 7, ...
                       'WhiteSpace', ' ', 'ConsecutiveDelimitersRule', 'join');
table4 = readtable('../OGLEIV/table4_central.txt',opts);

%----------------------------------
%Table 5. Detection efficiencies for the analyzed fields.
%-------------------------------------


VarNames_table5 = {'bin', 'log_tE', 'BLG500', 'BLG501', 'BLG504', 'BLG505', 'BLG506', 'BLG511', 'BLG512', 'BLG534', 'BLG611'};
VarTypes_table5 = {'double', 'double', 'double', 'double', 'double', 'double', 'double',  'double', 'double', 'double', 'double'}; 

opts = delimitedTextImportOptions('VariableNames',VarNames_table5,'VariableTypes',VarTypes_table5,...
                                'Delimiter',delimiter, 'DataLines', 6, ...
                       'WhiteSpace', ' ', 'ConsecutiveDelimitersRule', 'join');
eff_field = readtable('../OGLEIV/table5_central.txt',opts);

field_list = ["BLG500", "BLG501", "BLG504", "BLG505", "BLG506", "BLG511", "BLG512", "BLG534", "BLG611"];

%------------------------------------------------
%Table 5. Surface density of stars in OGLE-IV subfields calculated
%using image-level simulations.
%---------------------------------------------------
varNames = {'field', 'ra', 'dec', 'glon', 'glat', 'sigma18', 'sigma21', 'N18', 'N21'};
varTypes = {'char','double', 'double', 'double', 'double', ...
    'double', 'double', 'double', 'double'};

opts = delimitedTextImportOptions('VariableNames',varNames,'VariableTypes',varTypes,...
                                'Delimiter',delimiter, 'DataLines', 47, ...
                       'WhiteSpace', ' ', 'ConsecutiveDelimitersRule', 'join');
table5 = readtable('../OGLEIV/table5.dat',opts);


%-----------------------------------------------
%Table 7. Microlensing optical depth and event rates in the OGLE-IV
% fields (averaged over sources brighter than I=21).
%-----------------------------------------------

VarNames_table7 = {'field', 'glon', 'glat', 'tau', 'tau_err', 'gam', 'gam_err', 'gam_deg2',...
    'gam_deg2_err', 't_E_mean', 't_E_mean_err', 'N_events', 'N_stars'};
VarTypes_table7 = {'string', 'double', 'double', 'double', 'double', 'double', 'double', ...
    'double', 'double', 'double', 'double', 'double', 'double'}; 

opts = delimitedTextImportOptions('VariableNames',VarNames_table7,'VariableTypes',VarTypes_table7,...
                                'Delimiter',delimiter, 'DataLines', 25, ...
                       'WhiteSpace', ' ', 'ConsecutiveDelimitersRule', 'join');
table7 = readtable('../OGLEIV/table7.dat',opts);

%-----------------------------------------------
%Raw te
%-----------------------------------------------

VarNames_table_event = {'Name', 't_0', 't_E', 'u_0'};
VarTypes_table_event = {'string', 'double', 'double', 'double'}; 

opts = delimitedTextImportOptions('VariableNames',VarNames_table_event,'VariableTypes',VarTypes_table_event,...
                                'Delimiter',delimiter, 'DataLines', 8, ...
                       'WhiteSpace', ' ', 'ConsecutiveDelimitersRule', 'join');
table_event = readtable('../OGLEIV/OGLE-IV-events-FFP.txt',opts);

%-------------
% Choix du champ à analyser
%------------

% field = "BLG505";

for field = field_list
%-------------
%donnée expérimentales
%----------------

i_field = find(extractBetween(table_event.Name, 1,6) == field);
teff = table_event.t_E(i_field);

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

% on suit la direction du champ considéré

l = table6.glon(table6.field == field) *pi/180;    % direction d'observation en radian
b = table6.glat(table6.field == field) *pi/180;

uT = 1;		   % Seuil de d�tection en param�tre d'impact
AT = 3/sqrt(5);    % Seuil de d�tection en amplification


%-----------------------------------
% Param�tres de la fonction de masse
%-----------------------------------
global minf msup

minf=0.1;
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


%----------------------------------------
% recuperation des evenements selectionnes
%----------------------------------------
path_field = strcat('../graph/OGLEIV/', field);
mkdir(path_field);
mv evenements.txt path_field
mv simul_graph.txt path_field
end


%analyse des résultat
te_tot = [];
te_obs_tot =[];
te_obs_b_tot = [];

for field = field_list
    
path_field = strcat('../graph/OGLEIV/', field);
load(strcat(path_field, '/evenements.txt')

te = evenements(:,5);
te=te';

fid = fopen(stract(path_field,'/simul_para.txt');
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
% f = 0.05;   %fraction des évenements unblendé f = P(1)
% nbar = 4.51; % P(n) = fonction(nbar) = f avec P(n) la proba d'avoir n étoiles dans DeltaS

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

%---------------
% Données de l'expérience OGLE
%----------------------

%Exposure
exposure = 2011*table7.N_stars(table7.field == field) /365.25;
% exposure = 2741*sum(table5.N18(extractBetween(table5.field,1,6) == field))*1e-6 /365.25;

%Récupération des évènements du champs

disp(['exposure ogle = ', num2str(exposure)])
disp(['gamma observé ogle = ' num2str(table7.gam(find(extractBetween(table7.field, 1, 6) == field)))])
disp(['tau observé ogle = ' num2str(table7.tau(find(extractBetween(table7.field, 1, 6) == field)))])


%-------------
%Calcul de l'efficacité à partir des données de OGLE IV pour l'article de
%2017 (nature)
%-----------------
eff_2017 = eval(['eff_field.',convertStringsToChars(field)]);

te_inter = 10.^eff_field.log_tE;


teffmaxm=max(te_inter);
teffminm=min(te_inter);

i1_unblend = find((te<=teffmaxm)&(te>=teffminm));
i1_blend = find((teblend<=teffmaxm)&(teblend>=teffminm));

eff_unblend = zeros(1,length(te));	% applique une efficacite nulle aux durees superieures et inferieures
eff_blend = zeros(1,length(teblend));

eff_unblend(i1_unblend) = interp1(te_inter,eff_2017,te(i1_unblend));
eff_blend(i1_blend) = interp1(te_inter,eff_2017,teblend(i1_blend));


%---------------
%efficacité du dossier /eff
%------------------------
VarNames_eff_IV = {'log_tE_min', 'log_tE_max', 'efficiency'};
VarType_eff_IV = {'double', 'double', 'double'};

opts_eff = delimitedTextImportOptions('VariableNames',VarNames_eff_IV,'VariableTypes',VarType_eff_IV,...
                            'Delimiter',delimiter, 'DataLines', 5, ...
                   'WhiteSpace', ' ', 'ConsecutiveDelimitersRule', 'join');
               
eff_field_2019 = readtable(strcat("../OGLEIV/eff/", field, ".eff"),opts_eff);


% eff_2019 = eff_field_2019.efficiency;


te_inter_min = 10.^(eff_field_2019.log_tE_min);
te_inter_max = 10.^(eff_field_2019.log_tE_max);

teffmaxm=max(te_inter_max);
teffminm=min(te_inter_min);

% eff_unblend = zeros(1,length(te));	% applique une efficacite nulle aux durees superieures et inferieures
% eff_blend = zeros(1,length(teblend));
% 
% for i = 1:length(te_inter_min)
%     i1_unblend = find(te>=te_inter_min(i) & te<=te_inter_max(i));
%     i1_blend = find(teblend>=te_inter_min(i) & teblend<=te_inter_max(i));
%     
%     eff_unblend(i1_unblend) = ones(size(i1_unblend)) .* eff_field_2019.efficiency(i);
%     eff_blend(i1_blend) = ones(size(i1_blend)) .* eff_field_2019.efficiency(i);
%     
% end


%--------------------------------------------------------------------------------------------------------------------------
% compare le nombre aleatoire precedent a l'efficacite que l'on vient de calculer afin de decider si l'evt est garde ou non
%--------------------------------------------------------------------------------------------------------------------------

%choix de l'éfficacité
eff = eff_2017;

%tirage au sort pour l'efficacité

ra_unblend = rand(1,length(te))*max(eff);
ra_blend = rand(1,length(teblend))*max(eff);

% On choisit l'efficacité ici en prenant les bons indices i

i = find(ra_unblend-eff_unblend<=0); 
teobs = te(i);

ib = find(ra_blend-eff_blend<=0); 
teobsblend = teblend(ib); % On récupère les éléments qui sont soumis au blending avec le calcul d'avant


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

te_tot = [te_tot te];
te_obs_tot =[te_obs_tot teobs];
te_obs_b_tot = [te_obs_b_tot teobsblend];
end
%------------------------
% affichage des resultats
%------------------------

te = te_tot;
teobs = te_obs_tot;
teobsblend = te_obs_b_tot;
teff = table_event.t_E

%Paramètre graph
bin_max = 100;
nbre_bin = bin_max/2;

%Trace la distrib de te  pour le modèle et la courbe stockée localement
[hist, edges] = histcounts(te, nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');

%tracé distribution avec blending uniquement la courbe stockée localement
[histb, edges] = histcounts(teblend, nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');

%Courbe expérimentale (avec l'efficacité) :
[hist_obs, edges] = histcounts(teobs, nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');
[hist_obs_b, edges] = histcounts(teobsblend, nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');

centre = zeros(size(edges)-[0,1]);

for j =1:length(centre);
centre(j)=(edges(j)+edges(j+1))/2;
end

%expérience
[hist_exp_normalise, edges] = histcounts(teff, nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');
[hist_exp, edges] = histcounts(teff, nbre_bin, 'BinLimits',[0,bin_max]);

%graph en fonction de l'exposition
figure(17)
hold on;
plot(centre, hist_obs.*gamobs*exposure, 'red');
plot(centre, hist_obs_b*gamobsb*exposure, 'black');
plot(centre, hist_exp)
legend('hist modèle', 'hist modèle avec blending (f=0.5)', strcat('OGLE IV,  ', field))
xlabel('t_{e}')
ylabel('Nombre d''évènements par unité de t_{e}')

%Graph normalisé expérience et exp simulée avec l'efficacité pour OGLE IV
figure(18)
hold on;
plot(centre, hist_obs, 'red');
plot(centre, hist_obs_b, 'black');
plot(centre, hist_exp_normalise)
legend('hist modèle', 'hist modèle avec blending (f=0.5)', strcat('OGLE IV,  ', field))
xlabel('t_{e}')
ylabel('Nombre d''évènements par unité de t_{e}')


%%
%Graph normalisé expérience et exp simulée avec l'efficacité pour OGLE III
% figure(18)
% hold on;
% plot(centre, hist_obs, 'red');
% plot(centre, hist_obs_b, 'black');
% M = length(hist_exp_err);
% plot(edges(sort([1:M 1:M])), [0 , 0, hist_exp_err(sort([1:M 2:M-1]))])
% M = length(hist_exp_BW);
% plot(edges(sort([1:M 1:M])), [0 , 0, hist_exp_BW(sort([1:M 2:M-1]))], 'g')
% legend('hist modèle', 'hist modèle avec blending (f=0.5)', 'OGLE III (all stars)', 'OGLE III (fenêtre de Baade)')
% xlabel('t_{e}')
% ylabel('Nombre d''évènements par unité de t_{e}')

%graph noramlisé avec blending
% figure(1);
% hold on;
% plot(centre, hist, 'red')
% plot(centre, histb, 'black')
% title('Blending black et sans blending rouge)');
% hold off;
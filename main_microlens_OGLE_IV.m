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

%-----------------
%Table 3. Best-fitting parameters of the analyzed microlensing events in low-cadence OGLE fields.
%-----------------
varNames = {'name', 'field', 'star_id', 'ra', 'dec', 'ra_deg', 'dec_deg', 'glon',...
    'glat', 't0_best', 'tE_best', 'u0_best', 'Is_best', 'fs_best', 't0_med', 't0_err1', 't0_err2', 'tE_med',...
    'tE_err1', 'tE_err2', 'u0_med', 'u0_err1', 'u0_err2', 'Is_med', 'Is_err1', 'Is_err2', 'fs_med', 'fs_err1', 'fs_err2', 'weight', 'ews_id'} ;
varTypes = {'char','string', 'int32' ,'char', 'char', 'double', 'double', 'double', 'double', ...
    'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', ...
    'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'char'} ;

opts = delimitedTextImportOptions('VariableNames',varNames,'VariableTypes',varTypes,...
                                'Delimiter',delimiter, 'DataLines', 47, ...
                       'WhiteSpace', ' ', 'ConsecutiveDelimitersRule', 'join');
table3 = readtable('../OGLEIV/table3_corrected.dat',opts);

%-------------
% Choix du champ à analyser
%------------

% field = 'BLG535';
field_list = table7.field(find(table7.N_events>5));
% field_list = "BLG515";
for i_field = 55:56
%     :length(field_list)

field = field_list(i_field);
field = "BLG525";
disp(field)
% field = field_list(i_field);
disp(field)
global vlimit

vlimit = 1000e3;

%----------------------
% Nombre de simulations
%----------------------

n = 40e6;
nbsimul=1; %a augmenter pour meilleure stat

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
l = table6.glon(table6.field == field) *pi/180;    % direction d'observation en radian
b = table6.glat(table6.field == field) *pi/180;

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

%----------------------------------------
% recuperation des evenements selectionnes
%----------------------------------------
path_field = strcat('../graph/modif_FM/OGLEIV/', field);
movefile("evenements.txt", strcat(path_field, "/evenements_1.txt"))
movefile ("simul_para.txt", strcat(path_field, "/simul_para_1.txt"))
end


%%
%analyse des résultat
te_tot = [];
te_obs_tot =[];
te_obs_b_tot = [];
gam_obs_tot = [];
tau_obs_tot = [];
% field_list = ["BLG505"];

for i_field = 1:length(field_list)
field = field_list(i_field);
disp(field)

path_field = strcat('../graph/modif_FM/OGLEIV/', field);
evenements = load(strcat(path_field, '/evenements_2.txt'));
te = evenements(:,5);
te=te';

fid = fopen(strcat(path_field,'/simul_para_2.txt'));
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
% disp(['gamma (calcule par le te moyen) =    ' num2str(gamma)]);


gam1=4*sqrt(GMsol)/c*uT/sqrt(pc*pc*pc)*length(te)/(n*nbsimul)*86400*365.25*1e6/mmean;
gam=gam1*Gammax;
% disp(['gamma (integre par MC) = ' num2str(gam)]);


ttobs=tau/(gam/1e6/365.25);
% disp(['<tobs> (en jours) = ' num2str(ttobs)]);

taur=gam*pi/2*uT*mean(te)/365.25/1e6;
taur=real(taur);
% disp(['tau (avec gamma integré par MC) = ' num2str(taur)]);

%------------------------
% Application du blending
%------------------------

% Param�tre pour le blending 
f = 0.05;   %fraction des évenements unblendé f = P(1)
nbar = 4.51; % P(n) = fonction(nbar) = f avec P(n) la proba d'avoir n étoiles dans DeltaS

f = 0.5;
nbar = 1.257;
%retourn teblend (histogramme corrigé) et taurblend (profondeur optique corrigée)
script_blending 

%-----------
%Choix de l'expérience à analyser
%donne teff : tet des observations et eff : efficacité
%------------


temax = 100;
nbre_bin = temax;

exp_ogle_IV_2019


%---------------
%calcul de gamma
%---------------

% disp('Grandeurs avec intervention de l''efficacite experimentale :')

%------------------------------------------------------------------------------------------------
% on ne peut pas calculer le gamma par la formule avec le tobs, car le tau ne prend pas en compte
% l'efficacite. Par contre, on peut deduire tau experimental a partir du gamma calcule par MC
%------------------------------------------------------------------------------------------------

gamobs = gam/length(te)*length(teobs);
% disp(['gamma (integre par MC) = ' num2str(gamobs)]);
% 
tauobs=gamobs*pi/2*uT*mean(teobs)/365.25/1e6;
% disp(['tau obs (calcule par le te moyen) = ' num2str(tauobs)]);
% 
% tauobsblend=tauobs * gmean * (nbar/(1-exp(-nbar)));
% tauobsblend=real(tauobsblend);
% % disp(['tau observé avec blending (Alibert 2005)  = ' num2str(tauobsblend)]);
% 
gamobsb = gam/length(te)*length(teobsblend);
% disp(['gamma avec blending (integre par MC) = ' num2str(gamobsb)]);
% 
tauobsb=gamobsb*pi/2*uT*mean(teobsblend)/365.25/1e6;
% disp(['tau obs avec blending (calcule par le te moyen) = ' num2str(tauobsb)]);
% 
% disp(['rapport tau_blend/tau_obs_théorique = ' num2str(tauobsb/tauobs)]);
te_tot = [te_tot te];
te_obs_tot =[te_obs_tot teobs];
te_obs_b_tot = [te_obs_b_tot teobsblend];
gam_obs_tot = [gam_obs_tot, gamobs];
tau_obs_tot = [tau_obs_tot, tauobs];
end

%------------------------
% affichage des resultats
%------------------------

te = te_tot;
teobs = te_obs_tot;
teobsblend = te_obs_b_tot;

id_field = find(ismember(extractBetween(table3.field, 1, 6), field_list) & table3.Is_med<=21);
teff = table3.tE_best(id_field);

%prise en compte de l'efficacité pour l'histogramme corrigé
% N_events = table7.N_events(table7.field == field);
% gam_ogle = table7.gam(find(table7.field == field));
% mean_eff = N_events/(gam_ogle*exposure);

%------------------------
% affichage des resultats
%------------------------

%Paramètre graph
bin_max = 100;
nbre_bin = bin_max/2;


%Trace la distribde te  pour le modèle et la courbe stockée localement
[hist, edges] = histcounts(te, nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');

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

% %graph en fonction de l'exposition
% figure(17)
% hold on;
% plot(centre, hist_obs.*gamobs*exposure*mean_eff, 'red');
% plot(centre, hist_obs_b*gamobsb*exposure*mean_eff, 'black');
% histogram(teff, nbre_bin, 'BinLimits',[0,bin_max])
% legend('hist modèle', 'hist modèle avec blending (f=0.5)', strcat('OGLE IV,  ', field))
% xlabel('t_{e}')
% ylabel('Nombre d''évènements par unité de t_{e}')

%Graph normalisé expérience et exp simulée avec l'efficacité pour OGLE IV
if ishandle(18)
    close(18)
end
figure(18)
hold on;
plot(centre, hist_obs, 'red');
plot(centre, hist_obs_b, 'black');
M = length(hist_exp_err);
plot(edges(sort([1:M 1:M])), [0 , 0, hist_exp_err(sort([1:M 2:M-1]))])
legend('hist modèle', 'hist modèle avec blending (f=0.5)', 'OGLE IV')
xlabel('t_{e}')
ylabel('Nombre d''évènements par unité de t_{e}')

%Graph log
te_min = min(teff);
% te_min = 10^0.2;
te_max = max(teff); M = 26;
edges_log=te_min*(te_max/te_min).^([0:M]/M);
x=edges_log(sort([1:M+1 1:M+1])); 

%calcul centre pour errorbar
centre = zeros(size(edges_log)-[0,1]);
for j =1:length(centre);
centre(j)=(edges_log(j)+edges_log(j+1))/2;
end

hist_exp_log = histcounts(teff,edges_log);
hist_obs_log = histcounts(teobs,edges_log);
hist_obs_b_log = histcounts(teobsblend,edges_log);

if ishandle(2)
    close(2)
end
figure(2)
hold on
plot(centre, hist_obs_log./length(teobs))
plot(centre, hist_obs_b_log./length(teobsblend))
plot(x, [0 hist_exp_log(sort([1:M 1:M])) 0]./length(teff), 'b')
errorbar(centre, hist_exp_log./length(teff), 1./sqrt(hist_exp_log)./length(teff), 'b.')

% errorbar((edges_log + [edges_log(2:end) te_max+100])./2, hist_exp_log./length(teff), mean(teff)./sqrt(hist_exp_log)./length(teff), mean(teff)./sqrt(hist_exp_log)./length(teff), '.')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
legend('hist modèle', 'hist modèle avec blending (f=0.5)', 'OGLE IV')
legend('Location', 'best')
xlabel('log(t_{e})')
ylabel('Nombre d''évènements par unité de t_{e} (échelle log)')
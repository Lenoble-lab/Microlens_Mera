%But calculer tau, gamma et <te> pour chaque champ de OGLEIV selon un modèle considéré
%On charge les données OGLE IV et ensuite on calcule pour chaque direction (=champ) via MC ces données

clear

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


fichier ='../graph/modif_FM/comp_modele_OGLEIV.txt';

fid_model = fopen(fichier, 'w');

fprintf(fid_model, ['Calcul pour chaque champ des données de mon modèle, comparée avec celle de OGLE IV \n' ...
'Les variables stockées sont, dans l ordre : \n ' ...
'Champ; \n longitude; \n lattitude; \n N_stars; \n N_events; ' ...
'\n gam_OGLE; \n tau_OGLE (10^-6); \n mean_te_ogle; \n ' ...
'Gammax \n' ...
'gamma calculé par le modèle (brut); \n tau (10^-6); ----------------- \n mean_te ----------------------\n ' ...
'gamma avec efficacité du champ considéré \n tau (10^-6); ---------------- \n mean_te -------------------\n '...
'gamma avec efficacité et blending \n tau (10^-6); ---------------- \n mean_te --------------------\n' ]);


global vsr vsp vst vlp vlt vlr

%----------------------------------
% Constantes physiques (unites SI)
%----------------------------------
G=6.672e-11;	pc=3.08567802e16;
kpc=pc*1e3;  	Msol=1.989e30;
c=299792458;	GMsol=1.32712497e20;

%-------------------------------------------------
% Param�tres vlimit (vitesse perpendiculaire maxi)
%-------------------------------------------------

global vlimit

vlimit = 1000e3;

%----------------------
% Nombre de simulations
%----------------------

n = 10e5;
nbsimul=10; %a augmenter pour meilleure stat
nbMAX=500;

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

fprintf(fid_model, ' n = %10.0f, nbsimul = %10.0f, nbmax = %10.0f \n', [n, nbsimul, nbMAX]);


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


% for field = table7.field
for id_field = 1:length(table7.field)
 
field = table7.field(id_field);  
disp(field)
global l b

% on suit la direction du champ considéré

l = table7.glon(table7.field == field) *pi/180;    % direction d'observation en radian
b = table7.glat(table7.field == field) *pi/180;

uT = 1;		   % Seuil de d�tection en param�tre d'impact
AT = 3/sqrt(5);    % Seuil de d�tection en amplification


%-------------------
% Monte-Carlo
%-------------------
main

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
% f = 0.2;
% nbar = 2.6;
% f = 1;
% nbra = 0;
%retourn teblend (histogramme corrigé) et taurblend (profondeur optique corrigée)
script_blending 


%Calcul de l'efficacité
VarNames_eff_IV = {'log_tE_min', 'log_tE_max', 'efficiency'};
VarType_eff_IV = {'double', 'double', 'double'};

opts_eff = delimitedTextImportOptions('VariableNames',VarNames_eff_IV,'VariableTypes',VarType_eff_IV,...
                            'Delimiter',delimiter, 'DataLines', 5, ...
                   'WhiteSpace', ' ', 'ConsecutiveDelimitersRule', 'join');
               
eff_field = readtable(strcat("../OGLEIV/eff/", field, ".eff"),opts_eff);


eff = eff_field.efficiency;

te_inter_min = 10.^(eff_field.log_tE_min);
te_inter_max = 10.^(eff_field.log_tE_max);

teffmaxm=max(te_inter_max);
teffminm=min(te_inter_min);

eff_unblend = zeros(1,length(te));	% applique une efficacite nulle aux durees superieures et inferieures
eff_blend = zeros(1,length(teblend));

for i = 1:length(te_inter_min)
    i1_unblend = find(te>=te_inter_min(i) & te<=te_inter_max(i));
    i1_blend = find(teblend>=te_inter_min(i) & teblend<=te_inter_max(i));
    
    eff_unblend(i1_unblend) = ones(size(i1_unblend)) .* eff_field.efficiency(i);
    eff_blend(i1_blend) = ones(size(i1_blend)) .* eff_field.efficiency(i);
    
end

%tirage au sort pour l'efficacité
ra_unblend = rand(1,length(te))*max(eff_field.efficiency);
ra_blend = rand(1,length(teblend))*max(eff_field.efficiency);

% On choisit l'efficacité ici en prenant les bons indices i

i = find(ra_unblend-eff_unblend<=0); 
teobs = te(i);

ib = find(ra_blend-eff_blend<=0); 
teobsblend = teblend(ib); % On récupère les éléments qui sont soumis au blending avec le calcul d'avant

%-------------------------------
%Exploitation des résultats avec efficacité expérimentale
%-------------------------------

gamobs = gam/length(te)*length(teobs);
% disp(['gamma (integre par MC) = ' num2str(gamobs)]);

tauobs=gamobs*pi/2*uT*mean(teobs)/365.25/1e6;
% disp(['tau obs (calcule par le te moyen) = ' num2str(tauobs)]);

gamobsb = gam/length(te)*length(teobsblend);
% disp(['gamma avec blending (integre par MC) = ' num2str(gamobsb)]);

tauobsb=gamobsb*pi/2*uT*mean(teobsblend)/365.25/1e6;
% disp(['tau obs avec blending (calcule par le te moyen) = ' num2str(tauobsb)]);

%----------------------------------
%enregistrement dans le fichier txt
%------------------------------------

disp(gam)
evnts = [l*180/pi; b*180/pi; table7.N_stars(table7.field == field); table7.N_events(table7.field == field) ; ...
table7.gam(table7.field == field); table7.tau(table7.field == field); table7.t_E_mean(table7.field == field); ...
Gammax; gam; tau*1e6; mean(te); gamobs; tauobs*1e6; mean(teobs); ...
gamobsb; tauobsb*1e6; mean(teblend)];
fprintf(fid_model, field);
fprintf(fid_model,'  %12.8f  %12.8f  %12.8f  %12.8f %12.8f %12.8f  %12.8f  %12.8f %12.8f  %12.8f  %12.8f %12.8f  %12.8f  %12.8f %12.8f  %12.8f  %12.8f\n',evnts);
end

%--------------------------------------
%% 2nd partie : exploitation des résultats
%------------------------------

%-----------------------------------------------
%Table information about fields
%------------------------------------------------
delimiter = ' ';
VarNames_comp = {'field', 'glon', 'glat', 'N_stars', 'N_events', 'gam_OGLE', 'tau_OGLE', 'mean_te_ogle', 'Gammax', 'gamma_brut', 'tau_brut','mean_te_brut', 'gamma_eff', 'tau_eff', 'mean_te_eff', ...
    'gamma_eff_blend', 'tau_eff_blend', 'mean_te_eff_blend'};


VarTypes_comp = {'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}; 

opts = delimitedTextImportOptions('VariableNames',VarNames_comp,'VariableTypes',VarTypes_comp,...
                                'Delimiter',delimiter, 'DataLines', 22, ...
                       'WhiteSpace', ' ', 'ConsecutiveDelimitersRule', 'join');
comp_modele = readtable('../graph/comp_modele_OGLEIV.txt',opts);


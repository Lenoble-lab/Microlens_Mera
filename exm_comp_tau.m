%But tracé tau en fonction de b pour le modèle et pour deifférentes
%expériences
clear; close all;
%-----------------------------------------------
%OGLE IV
%Table 7. Microlensing optical depth and event rates in the OGLE-IV
% fields (averaged over sources brighter than I=21).
%-----------------------------------------------
delimiter = ' ';
VarNames_table7 = {'field', 'glon', 'glat', 'tau', 'tau_err', 'gam', 'gam_err', 'gam_deg2',...
    'gam_deg2_err', 't_E_mean', 't_E_mean_err', 'N_events', 'N_stars'};
VarTypes_table7 = {'string', 'double', 'double', 'double', 'double', 'double', 'double', ...
    'double', 'double', 'double', 'double', 'double', 'double'}; 

opts = delimitedTextImportOptions('VariableNames',VarNames_table7,'VariableTypes',VarTypes_table7,...
                                'Delimiter',delimiter, 'DataLines', 25, ...
                       'WhiteSpace', ' ', 'ConsecutiveDelimitersRule', 'join');
table7 = readtable('../OGLEIV/table7.dat',opts);

%-----------------------------------------------
% MOA 2016
% Table 4. Average microlensing optical depth and event rates at the position 
%of each subfield for the all-source sample
%------------------------------------------------

delimiter = ' ';
VarNames_table4 = {'Field', 'glon', 'glat',  'Nsub',  'Nstar',  'Nev', 'tau',  'etau', 'e_tau', 'E_tau', 'Gamma', 'eGamma', ...
'e_Gamma', 'E_Gamma', 'Gammad', 'eGammad', 'e_Gammad', 'E_Gammad'};
VarTypes_table4 = {'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double',...
 'double', 'double', 'double', 'double', 'double', 'double', 'double'}; 

opts = delimitedTextImportOptions('VariableNames',VarNames_table4,'VariableTypes',VarTypes_table4,...
                                'Delimiter',delimiter, 'DataLines', 37, ...
                       'WhiteSpace', ' ', 'ConsecutiveDelimitersRule', 'join');
table4 = readtable('../sumi_penny/table4.txt',opts);

%-----------------------------------------------
%Data from my model on OGLE IV field
%------------------------------------------------
delimiter = ' ';
VarNames_comp = {'field', 'glon', 'glat', 'N_stars', 'N_events', 'gam_OGLE', 'tau_OGLE', 'mean_te_ogle', 'Gammax', ...
    'gamma_brut', 'tau_brut','mean_te_brut', 'gamma_eff', 'tau_eff', 'mean_te_eff', ...
    'gamma_eff_blend', 'tau_eff_blend', 'mean_te_eff_blend'};


VarTypes_comp = {'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}; 

opts = delimitedTextImportOptions('VariableNames',VarNames_comp,'VariableTypes',VarTypes_comp,...
                                'Delimiter',delimiter, 'DataLines', 22, ...
                       'WhiteSpace', ' ', 'ConsecutiveDelimitersRule', 'join');
comp_modele = readtable('../graph/modif_FM/comp_modele_OGLEIV.txt',opts);

%----------------
%Tracé profondeur optique en fonction de la lattitude galactique au centre
%---------------------

glong = 0;
bin_long = 3;
bin_lat = 0.1;

%MOA II
glat = (0:bin_lat:1) .*-6 -1;

tau_mean_MOA = zeros(size(glat)-1);
lat_mean_MOA = zeros(size(glat)-1);
e_tau_mean_MOA = zeros(size(glat)-1);
E_tau_mean_MOA = zeros(size(glat)-1);

for i = 2:length(glat)
i0 = find(abs(table4.glon - glong)<bin_long & glat(i-1)>table4.glat & table4.glat>glat(i));
i = i-1;
tau_mean_MOA(i) = mean(table4.tau(i0));
lat_mean_MOA(i) = mean(table4.glat(i0));
e_tau_mean_MOA(i) = mean(table4.e_tau(i0));
E_tau_mean_MOA(i) = mean(table4.E_tau(i0));
end


%ogle IV et comp modèle
i0 = find(abs(table7.glon - glong)<bin_long & table7.glat<-1);

lat_mean_IV = zeros(size(glat)-1);


tau_mean_IV = zeros(size(glat)-1);
e_tau_mean_IV = zeros(size(glat)-1);
tau_brut = zeros(size(glat)-1);
tau_eff = zeros(size(glat)-1);
tau_eff_blend = zeros(size(glat)-1);

gam_mean_IV = zeros(size(glat)-1);
e_gam_mean_IV = zeros(size(glat)-1);
gam_brut = zeros(size(glat)-1);
gam_eff = zeros(size(glat)-1);
gam_eff_blend = zeros(size(glat)-1);

te_mean_IV = zeros(size(glat)-1);
e_te_mean_IV = zeros(size(glat)-1);
te_brut = zeros(size(glat)-1);
te_eff = zeros(size(glat)-1);
te_eff_blend = zeros(size(glat)-1);


for i = 2:length(glat)
i0 = find(abs(table7.glon - glong)<bin_long & glat(i-1)>table7.glat & table7.glat>glat(i));
i = i-1;
tau_mean_IV(i) = mean(table7.tau(i0));
lat_mean_IV(i) = mean(table7.glat(i0));
e_tau_mean_IV(i) = mean(table7.tau_err(i0));
gam_mean_IV(i) = mean(table7.gam(i0));
e_gam_mean_IV(i) = mean(table7.gam_err(i0));
te_mean_IV(i) = mean(table7.t_E_mean(i0));
e_te_mean_IV(i) = mean(table7.t_E_mean_err(i0));

tau_brut(i) = mean(comp_modele.tau_brut(i0));
tau_eff(i) = mean(comp_modele.tau_eff(i0));
tau_eff_blend(i) = mean(comp_modele.tau_eff_blend(i0));

gam_brut(i) = mean(comp_modele.gamma_brut(i0));
gam_eff(i) = mean(comp_modele.gamma_eff(i0));
gam_eff_blend(i) = mean(comp_modele.gamma_eff_blend(i0));

te_brut(i) = mean(comp_modele.mean_te_brut(i0));
te_eff(i) = mean(comp_modele.mean_te_eff(i0));
te_eff_blend(i) = mean(comp_modele.mean_te_eff_blend(i0));
end


figure(1)
hold on
errorbar(lat_mean_MOA, tau_mean_MOA, e_tau_mean_MOA, E_tau_mean_MOA, 'o')
errorbar(lat_mean_IV, tau_mean_IV, e_tau_mean_IV, 'o')
plot(lat_mean_IV, tau_brut)
plot(lat_mean_IV, tau_eff)
plot(lat_mean_IV, tau_eff_blend)

legend('MOA II', 'OGLE IV', 'Modèle brut', 'Modèle avec efficacité expérimentale (OGLE IV)', 'Modèle avec efficacité expérimentale et blending (OGLE IV)')
xlabel('b (deg)')
ylabel('\tau \times 10^{-6}')

figure(2)
hold on
errorbar(lat_mean_IV, gam_mean_IV, e_gam_mean_IV, 'o')
plot(lat_mean_IV, gam_brut)
plot(lat_mean_IV, gam_eff)
plot(lat_mean_IV, gam_eff_blend)

legend('OGLE IV', 'Modèle brut', 'Modèle avec efficacité expérimentale (OGLE IV)', 'Modèle avec efficacité expérimentale et blending (OGLE IV)')
xlabel('b (deg)')
ylabel('\Gamma')

figure(3)
hold on
errorbar(lat_mean_IV, te_mean_IV, e_te_mean_IV, 'o')
plot(lat_mean_IV, te_brut)
plot(lat_mean_IV, te_eff)
plot(lat_mean_IV, te_eff_blend)

legend('OGLE IV', 'Modèle brut', 'Modèle avec efficacité expérimentale (OGLE IV)', 'Modèle avec efficacité expérimentale et blending (OGLE IV)')
xlabel('b (deg)')
ylabel('<t_{e}>')
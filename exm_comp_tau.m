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


%----------------
%Tracé profondeur optique en fonction de la lattitude galactique au centre
%---------------------

glong = 0;
bin_long = 3;
bin_lat = 0.1;

%MOA II
glat = (0:bin_lat:1) .*-6 -1;

tau_mean = zeros(size(glat)-1);
lat_mean = zeros(size(glat)-1);
e_tau_mean = zeros(size(glat)-1);
E_tau_mean = zeros(size(glat)-1);

for i = 2:length(glat)
i0 = find(abs(table4.glon - glong)<bin_long & glat(i-1)>table4.glat & table4.glat>glat(i));
i = i-1;
tau_mean(i) = mean(table4.tau(i0));
lat_mean(i) = mean(table4.glat(i0));
e_tau_mean(i) = mean(table4.e_tau(i0));
E_tau_mean(i) = mean(table4.E_tau(i0));
end


%ogle IV
i0 = find(abs(table7.glon - glong)<bin_long & table7.glat<-1);

tau_mean_IV = zeros(size(glat)-1);
lat_mean_IV = zeros(size(glat)-1);
e_tau_mean_IV = zeros(size(glat)-1);

for i = 2:length(glat)
i0 = find(abs(table7.glon - glong)<bin_long & glat(i-1)>table7.glat & table7.glat>glat(i));
i = i-1;
tau_mean_IV(i) = mean(table7.tau(i0));
lat_mean_IV(i) = mean(table7.glat(i0));
e_tau_mean_IV(i) = mean(table7.tau_err(i0));
end

%Calcul modèle
tau_load = load('graph_iso_model.mat');

[L, B] = meshgrid(tau_load.L, tau_load.B);

i1 = find(abs(tau_load.L)<2);
i2 = find(tau_load.B<-1 & tau_load.B>-8);
i_long = find(abs(tau_load.L-glong) == 0);

figure(1)
hold on
errorbar(lat_mean, tau_mean, e_tau_mean, E_tau_mean, 'o')
errorbar(lat_mean_IV, tau_mean_IV, e_tau_mean_IV, 'o')
plot(tau_load.B(i2), mean(tau_load.tau_table(i2,i1),2)*1e6)
% plot(tau_load.B(i2), tau_load.tau_table(i2,i_long)*1e6)

legend('MOA II', 'OGLE IV', 'Modèle')
xlabel('b (deg)')
ylabel('\tau \times 10^{-6}')

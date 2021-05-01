% Prgm pour tracer les obs de Mroz 2020 ds le plan galactique
% Table 2 : tau, Gamma et <te> en fonction longitude galactique

%-----------------------------------------------
%Table 2. Microlensing Optical Depth and Event Rate toward the Galactic Plane
%------------------------------------------------
delimiter = ' ';
VarNames_table6 = {'l_min', 'l_max', 'b_min', 'b_max', 'nevents', 'nstars', 'tau', 'e_tau', 'gamma', 'e_gamma', 'gamma_d', 'e_gamma_d', 'te_mean', 'e_te_mean'};
VarTypes_table6 = {'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}; 

opts = delimitedTextImportOptions('VariableNames',VarNames_table6,'VariableTypes',VarTypes_table6,...
                                'Delimiter',delimiter, 'DataLines', 26, ...
                       'WhiteSpace', ' ', 'ConsecutiveDelimitersRule', 'join');
table2 = readtable('../Mroz_2020/table2.txt',opts);

%--------------------------
% load result local model
%---------------------------
local_model = load('graph_iso_model.mat');

% Convertion radian->degré
local_model.B_table = local_model.B_table * 180/pi;
local_model.L_table = local_model.L_table * 180/pi;

local_tau = mean(local_model.tau_table, 2);
local_gamma = mean(local_model.gamma_table, 2);

centre = (table2.l_min+table2.l_max)/2;
centre(1:8) = centre(1:8)-360;

close all;

%tau comparaison longitude
figure(1)
hold on

%plot(centre(1:11), table2.tau(1:11), 'ob')
errorbar(centre(1:11), table2.tau(1:11), table2.e_tau(1:11))
plot(local_model.L_table, local_tau*1e6)
legend('obs OGLE IV', 'modèle')
grid on;
grid minor
xlabel('longitude galactique (en degrès)')
ylabel('tau (x1e-6)')
axis([-200 200 0 0.8])
%gamma
figure(2)
hold on

%plot(centre(1:11), table2.gamma(1:11), 'b')
errorbar(centre(1:11), table2.gamma(1:11), table2.e_gamma(1:11))
plot(local_model.L_table, local_gamma)

axis([-200 200 0 4])
grid on;
grid minor
legend('obs OGLE IV', 'modèle')
xlabel('longitude galactique (en degrès)')
ylabel('gamma')

%% Lattitude

L_table = mod(local_model.L_table, 360);

i_lat = find(L_table<330 & L_table>240);

tau_lat = mean(local_model.tau_table(i_lat,:), 1);
gamma_lat = mean(local_model.gamma_table(i_lat,:), 1);

centre = (table2.b_min(12:18)+table2.b_max(12:18))/2;

close all;

%tau comparaison longitude
figure(1)
hold on

%plot(centre(1:11), table2.tau(1:11), 'ob')
errorbar(centre, table2.tau(12:18), table2.e_tau(12:18))
plot(local_model.B_table, tau_lat*1e6)
legend('obs OGLE IV', 'modèle')
grid on;
grid minor
xlabel('latitude galactique (en degrès)')
ylabel('tau (x1e-6)')
% axis([-200 200 0 0.8])
%gamma
figure(2)
hold on

%plot(centre(1:11), table2.gamma(1:11), 'b')
errorbar(centre, table2.gamma(12:18), table2.e_gamma(12:18))
plot(local_model.B_table, gamma_lat)

% axis([-200 200 0 4])
grid on;
grid minor
legend('obs OGLE IV', 'modèle')
xlabel('latitude galactique (en degrès)')
ylabel('gamma')
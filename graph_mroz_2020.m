% Prgm pour tracer les obs de Mroz 2020 ds le plan galactique
% Table 2 : tau, Gamma et <te> en fonction longitude galactique

%-----------------------------------------------
%Table 2. Microlensing Optical Depth and Event Rate toward the Galactic Plane
%------------------------------------------------
delimiter = ' ';
VarNames_table6 = {'l_min', 'l_max', 'b_min', 'b_max', 'nevents', 'nstars', 'tau', 'e_tau', 'gamma', 'e_gamma', 'gamma_d', 'e_gamma_d', 'te_mean', 'e_te_mean'};
VarTypes_table6 = {'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}; 

opts = delimitedTextImportOptions('VariableNames',VarNames_table6,'VariableTypes',VarTypes_table6,...
                                'Delimiter',delimiter, 'DataLines', 26, ...
                       'WhiteSpace', ' ', 'ConsecutiveDelimitersRule', 'join');
table2 = readtable('../meroz_2020/table2.txt',opts);

%--------------------------
% load result local model
%---------------------------
local_model = load('graph_iso_model.mat');

% Convertion radian->degré
local_model.B_table = tau_load.B_table * 180/pi;
local_model.L_table = tau_load.L_table * 180/pi;

local_tau = mean(local_model.tau_table, 2)
local_gamma = mean(local_model.tau_gamma, 2)

centre = (table2.l_min+table2.l_max)/2

%tau
figure(1)
hold on

plot(centre, table2.tau, 'b')
errorbar(centre, table2.tau, table2.e_tau, 'b.')
plot(local_model.L_table, local_tau*1e6)

grid on;
grid minor
xlabel('longitude galactique (en degrès)')
ylabel('tau (x1e-6)')

%gamma
figure(2)
hold on

plot(centre, table2.gamma, 'b')
errorbar(centre, table2.gamma, table2.e_gamma, 'b.')
plot(local_model.L_table, local_gamma)

grid on;
grid minor

xlabel('longitude galactique (en degrès)')
ylabel('tau (x1e-6)')
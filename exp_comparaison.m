%Tracé de courbe avec différente expériences
clear; close all
%-----------------
%Comparaison histogramme OGLE III et OGLE IV
%----------------

%---------------------
%Load table 2 de OGLE III
%----------------------

delimiter = ' ';
VarNames_table2 = {'id', 'ra', 'dec', 'field', 'N_stars', 't0', 'e_t0', 'E_t0', 'te', 'e_te', 'E_te', ...
'u0', 'e_u0', 'E_u0', 'fs', 'e_fs', 'E_fs', 'I0', 'e_I0', 'E_I0', 'Chi', 'Ndof','V' ,'e_V', 'EWS_id', 'glon', 'glat'};
VarTypes_table2 = {'double', 'string', 'string', 'string', 'double', 'double', 'double', ...
'double', 'double', 'double', 'double','double', 'double', 'double', 'double', ...
'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'string', 'double', 'double'}; 

opts = delimitedTextImportOptions('VariableNames',VarNames_table2,'VariableTypes',VarTypes_table2,...
                                'Delimiter',delimiter, 'DataLines', 8, ...
                       'WhiteSpace', ' ', 'ConsecutiveDelimitersRule', 'join');
table2_III = readtable('../OGLEIII/table2.txt',opts);

%-----------------
%Table 3 de OGLE IV
%-----------------

varNames = {'name', 'field', 'star_id', 'ra', 'dec', 'ra', 'dec', 'glon',...
    'glat', 't0_best', 'tE_best', 'u0_best', 'Is_best', 'fs_best', 't0_med', 't0_err1', 't0_err2', 'tE_med',...
    'tE_err1', 'tE_err2', 'u0_med', 'u0_err1', 'u0_err2', 'Is_med', 'Is_err1', 'Is_err2', 'fs_med', 'fs_err1', 'fs_err2', 'weight', 'ews_id'} ;
varTypes = {'char','string', 'int32' ,'char', 'char', 'double', 'double', 'double', 'double', ...
    'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', ...
    'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'char'} ;

opts = delimitedTextImportOptions('VariableNames',varNames,'VariableTypes',varTypes,...
                                'Delimiter',delimiter, 'DataLines', 47, ...
                       'WhiteSpace', ' ', 'ConsecutiveDelimitersRule', 'join');
table3_IV = readtable('../OGLEIV/table3.dat',opts);

%-----------------------
%Table 4 Number of events detected in individual timescale bins, OGLE IV
%2017
%------------------------------

delimiter = ' ';
VarNames_table4 = {'bin', 'log_tE', 'BLG500', 'BLG501', 'BLG504', 'BLG505', 'BLG506', 'BLG511', 'BLG512', 'BLG534', 'BLG611'};
VarTypes_table4 = {'double', 'double', 'double', 'double', 'double', 'double', 'double',  'double', 'double', 'double', 'double'}; 

opts = delimitedTextImportOptions('VariableNames',VarNames_table4,'VariableTypes',VarTypes_table4,...
                                'Delimiter',delimiter, 'DataLines', 7, ...
                       'WhiteSpace', ' ', 'ConsecutiveDelimitersRule', 'join');
table4_2017 = readtable('../OGLEIV/table4_central.txt',opts);

%---------------------------------------------------
%Table 4 : Microlensing events used in the optical depth and event rate measurements.
%-----------------------------------------------------
VarNames_table4_MOA = {'blank',  'ID', 'RA', 'Dec', 'Ndata', 't0', 'tE', 'e_tE', 'u0', 'e_u0', 'Is', 'chi2dof'};
VarTypes_table4_MOA = {'string', 'string', 'string', 'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'};

opts = delimitedTextImportOptions('VariableNames',VarNames_table4_MOA,'VariableTypes',VarTypes_table4_MOA,...
                                'Delimiter',delimiter, 'DataLines', 22, ...
                       'WhiteSpace', '  ', 'ConsecutiveDelimitersRule', 'join');
                   
table4_MOA = readtable('../MOA_II/Table4.dat',opts);


%--------------------------
%Efficacité
%---------------------------

%OGLE III
VarNames_eff = {'id', 't_e', 'N_events', 'undef', 'efficiency'};
VarTypes_eff = {'double', 'double', 'double', 'double', 'double'};

opts = delimitedTextImportOptions('VariableNames',VarNames_eff,'VariableTypes',VarTypes_eff,...
                                'Delimiter',delimiter,...
                       'WhiteSpace', ' ', 'ConsecutiveDelimitersRule', 'join');
eff_III = readtable('../OGLEIII/efficiency.dat',opts);


%OGLE IV

fields = ["BLG501", "BLG506", "BLG513", "BLG524"];

VarNames_eff_IV = {'log_tE_min', 'log_tE_max', 'efficiency'};
VarType_eff_IV = {'double', 'double', 'double'};


% opts_eff = delimitedTextImportOptions('VariableNames',VarNames_eff_IV,'VariableTypes',VarType_eff_IV,...
%         'Delimiter',' ', 'Data",Lines', 4, 'WhiteSpace', ' ', 'ConsecutiveDelimitersRule', 'join');
opts_eff = delimitedTextImportOptions('VariableNames',VarNames_eff_IV,'VariableTypes',VarType_eff_IV,...
                            'Delimiter',delimiter, 'DataLines', 5, ...
                   'WhiteSpace', ' ', 'ConsecutiveDelimitersRule', 'join');

% eff = readtable(strcat('../OGLEIV/eff/',extractBetween(table3.field(1), 1, 6),'.eff'),opts);

for i =1:length(fields)
assignin('base', ['eff_' num2str(i)], readtable(strcat("../OGLEIV/eff/", fields(i), ".eff"),opts_eff))
end

%Plot des graphes

%Comparaison différents champs de OGLE IV

field_2017 = ["BLG501", "BLG506"];
field_2019 = ["BLG535", "BLG513"];
%les bins sont imposés par le choix de OGLE IV
te_max = 150;

i_inf = find(10.^(table4_2017.log_tE+0.07) < te_max);
edges = 10.^[-1; table4_2017.log_tE(i_inf)+0.07];
i_ogle = find(table4_2017.log_tE<log10(max(edges)));

figure(1)
hold on;
for field = field_2017
hist_exp = eval(['table4_2017.', convertStringsToChars(field)]);
M = length(i_ogle);
% plot(10.^[-2;table4_2017.log_tE(i_ogle(sort([1:M 1:M])))], [0;0; hist_exp(sort([1:M-1 1:M-1])) ; 0], '-')
plot(10.^table4_2017.log_tE(i_ogle), hist_exp(i_ogle))
end

centre = zeros(size(edges)-[1,0]);
for j =1:length(centre);
centre(j)=(edges(j)+edges(j+1))/2;
end

for field = field_2019
id_field = find(extractBetween(table3_IV.field, 1, 6) == field & table3_IV.Is_med<=21);
teff = table3_IV.tE_best(id_field);
hist = histcounts(teff, edges);
M = length(hist);
% plot(edges(sort([1:M 1:M])), [0 , 0, hist(sort([1:M 2:M-1]))])
plot(centre, hist)
end
legend([field_2017, field_2019])
    

%%
%Histogramme complet des évènements
nbre_bin = 30;
figure(1)
semilogxhist(table3_IV.tE_best, nbre_bin)

hold on
semilogxhist(table2_III.te, nbre_bin)
semilogxhist(table4_MOA.tE, nbre_bin)

legend('OGLE IV', 'OGLE III', 'MOA')
xlabel('t_{e}')
ylabel('Nombre d''évènements par unité de t_{e}')


%Comparaison efficacité
figure(2)
%OGLE III
M = length(eff_III.t_e);
loglog(eff_III.t_e(sort([1:M 1:M])), [0 ; eff_III.efficiency(sort([1:M-1 1:M-1])) ; 0], '-')
hold on

%OGLE IV
for i =1:length(fields)
eff = eval(['eff_', num2str(i)]);
M = length(eff.log_tE_min);
te_graph = [eff.log_tE_min(1) ; eff.log_tE_min ; eff.log_tE_max; eff.log_tE_max(M)];
eff_graph = [0 ; eff.efficiency(sort([1:M 1:M])) ; 0];
loglog(10.^sort(te_graph), eff_graph)
end    

legend(['OGLE III', fields])
xlabel('t_{e}')
ylabel('Efficacité de détection')




function semilogxhist(val,M)
% semilogxhist - generate histogram with M bars and log-scale x axis
vmin=min(val); vmax=max(val);
edges=vmin*(vmax/vmin).^([0:M]/M);
count=histc(val,edges); 
if size(count,2)==1, count=count'; end 
x=edges(sort([1:M 1:M])); 
y=[0 count(sort([1:M-1 1:M-1])) 0];

% outline only: semilogx(x, y, '-');
plot(x, y, '-'); 
% fill(x, y, 'b'); 
set(gca,'XScale','log');
end
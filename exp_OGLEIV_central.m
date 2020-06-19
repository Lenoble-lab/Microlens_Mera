%Analyse des 9 champs centraux d'OGLE IV, ceux de l'article de Nature
%Champs observés très fréquement
%Mroz, P., Udalski, A., Skowron, J., et al. 2017, Nature 548, 183.

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
eff = readtable('../OGLEIV/table5_central.txt',opts);

field_list = ["BLG500", "BLG501", "BLG504", "BLG505", "BLG506", "BLG511", "BLG512", "BLG534", "BLG611"];
%-------------------
%Histogram of the diff fields
%------------------
nbre_bin = 30;
figure(1)
hold on
for i =1:9
plot(10.^table4.log_tE, eval(['table4.',convertStringsToChars(field_list(i))]));
end
set(gca, 'XScale', 'log')
legend(field_list)
xlabel('t_{e}')
ylabel('Nombre d''évènements par unité de t_{e}')

%------------------
%Efficiency
%--------------------------
figure(2)
hold on
for i =1:9
plot(10.^table4.log_tE, eval(['eff.',convertStringsToChars(field_list(i))]));
end
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
legend(field_list)
xlabel('t_{e}')
ylabel('Nombre d''évènements par unité de t_{e}')
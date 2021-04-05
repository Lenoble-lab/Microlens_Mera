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

%-----------------------------------------------
%Raw te
%-----------------------------------------------

VarNames_table_event = {'Name', 't_0', 't_E', 'u_0'};
VarTypes_table_event = {'string', 'double', 'double', 'double'}; 

opts = delimitedTextImportOptions('VariableNames',VarNames_table_event,'VariableTypes',VarTypes_table_event,...
                                'Delimiter',delimiter, 'DataLines', 8, ...
                       'WhiteSpace', ' ', 'ConsecutiveDelimitersRule', 'join');
table_event = readtable('../OGLEIV/OGLE-IV-events-FFP.txt',opts);

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


field_list = ["BLG500", "BLG501", "BLG504", "BLG505", "BLG506", "BLG511", "BLG512", "BLG534", "BLG611"];
%-------------------
%Histogram of the diff fields
%------------------
close all;

%Graph log
te_min = min(teff); 
% te_min = 0.1;
te_max = max(teff); M = 26;
edges_log=te_min*(te_max/te_min).^([0:M]/M);
x=edges_log(sort([1:M+1 1:M+1])); 

%calcul centre pour errorbar
centre = zeros(size(edges_log)-[0,1]);
for j =1:length(centre);
centre(j)=(edges_log(j)+edges_log(j+1))/2;
end 
legend_graph = [];
figure(1)
hold on
for i =1:2:9
teff = table_event.t_E(find(extractBetween(table_event.Name, 1, 6) == field_list(i)));

hist_exp_log = histcounts(teff,edges_log);
% plot(x, [0 hist_exp_log(sort([1:M 1:M])) 0]./length(teff))
plot(centre, hist_exp_log./length(teff))

% errorbar(centre, hist_exp_log./length(teff), 1./sqrt(hist_exp_log)./length(teff))
legend_graph = [legend_graph, field_list(i)];
end
plot(centre, histcounts(table3.tE_best, edges_log)/length(table3.tE_best), '--', 'LineWidth',2);
set(gca, 'XScale', 'log')
legend([legend_graph, 'OGLE IV total'])
legend('Location', 'best')
xlabel('log(t_{e})')
ylabel('Nombre d''évènements par unité de t_{e} (échelle log)')

%%
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
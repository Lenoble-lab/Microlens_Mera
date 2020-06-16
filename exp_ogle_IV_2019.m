%---------------------------------------------------------------------------------------------------------------
%Microlensing optical depth and event rate in the OGLE-IV Galactic plane fields
%
% P. Mroz, A. Udalski, M.K. Szymanski, I. Soszynski, P. Pietrukowicz,
% S. Kozlowski, J. Skowron, R. Poleski, K. Ulaczyk, M. Gromadzki,
% K. Rybicki, P. Iwanek, and M. Wrona
%---------------------------------------------------------------------------------------------------------------
% close all
% clear


% ttobs=tau/(gam/1e6/365.25);
% disp(['<tobs> (en jours) = ' num2str(ttobs)]);
% N=gam*exposure;
% disp(['nb d''evt  = ' num2str(N)]);

% taur=gam*pi/2*uT*mean(te)/365.25/1e6;
% taur=real(taur);
% disp(['tau (avec gamma integré par MC) = ' num2str(taur)]);

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
table3 = readtable('../OGLEIV/table3.dat',opts);


% centre = zeros(size(edges)-[0,1]);
% for j =1:length(centre);
% centre(j)=(edges(j)+edges(j+1))/2;
% end

% figure(1)
% loglog(centre, hist_teff)
% range = 0:1e-1*2:1 * log(400);
% bin_range = exp(range);
% 
% [bincounts] = histc(teff, bin_range);
% figure(2)
% loglog(bin_range, bincounts)

%----------------
%Tracé profondeur optique en fonction de la lattitude galactique au centre
%---------------------
%Longitude considérée
% long = 0;
% 
% i0 = find(abs(table6.glon - long)<0.5 & table6.glat<-1);
% 
% %Calcul modèle
% tau_load = load('graph_iso_model.mat');
% 
% [L, B] = meshgrid(tau_load.L, tau_load.B);
% 
% i1 = find(abs(tau_load.L)<2);
% i2 = find(tau_load.B<-1);
% i_long = find(abs(tau_load.L-long) == 0);
% 
% figure(1)
% hold on
% errorbar(table6.glat(i0), table7.tau(i0), table7.tau_err(i0), 'o')
% plot(tau_load.B(i2), mean(tau_load.tau_table(i2,i1),2)*1e6)
% plot(tau_load.B(i2), tau_load.tau_table(i2,i_long)*1e6)
% 
% legend('Mesure d''OGLE IV', 'Modèle')
% xlabel('b (deg)')
% ylabel('\tau \times 10^{-6}')

%---------------
%Efficacité
%-----------
% VarNames = {'log_tE_min', 'log_tE_max', 'efficiency'};
% VarType = {'double', 'double', 'double'};
% 
% opts = delimitedTextImportOptions('VariableNames',VarNames,'VariableTypes',VarType,...
%         'Delimiter',' ', 'DataLines', 5, 'WhiteSpace', ' ', 'ConsecutiveDelimitersRule', 'join');
%     
% % eff = readtable(strcat('../OGLEIV/eff/',extractBetween(table3.field(1), 1, 6),'.eff'),opts);
% eff = readtable(strcat('../OGLEIV/eff/','BLG500','.eff'),opts);

%----------------------------------------------
%Comparaison avec calculs du MC
%-----------------------------------------------

%Choix d'un champ 
field = "BLG612";

%Exposure
exposure = 2741*table7.N_stars(find(table7.field == field)) /365.25;
%Récupération des évènements du champs

disp(['exposure ogle = ', num2str(exposure)])
disp(['gamma observé ogle = ' num2str(table7.gam(find(extractBetween(table7.field, 1, 6) == field)))])
disp(['tau observé ogle = ' num2str(table7.tau(find(extractBetween(table7.field, 1, 6) == field)))])

id_field = find(extractBetween(table3.field, 1, 6) == field & table3.Is_med<=21);
teff = table3.tE_best(id_field);

%Calcul de l'efficacité
VarNames_eff_IV = {'log_tE_min', 'log_tE_max', 'efficiency'};
VarType_eff_IV = {'double', 'double', 'double'};

opts_eff = delimitedTextImportOptions('VariableNames',VarNames_eff_IV,'VariableTypes',VarType_eff_IV,...
                            'Delimiter',delimiter, 'DataLines', 5, ...
                   'WhiteSpace', ' ', 'ConsecutiveDelimitersRule', 'join');
               
eff_field = readtable(strcat("../OGLEIV/eff/", field, ".eff"),opts_eff);


eff = eff_field.efficiency;

figure(1)
loglog(table3.tE_best(id_field), table3.weight(id_field).^-1, 'o')
hold on

M = length(eff_field.log_tE_min);
te_graph = [eff_field.log_tE_min(1) ; eff_field.log_tE_min ; eff_field.log_tE_max; eff_field.log_tE_max(M)];
eff_graph = [0 ; eff_field.efficiency(sort([1:M 1:M])) ; 0];
loglog(10.^sort(te_graph), eff_graph)


%-----------------------------------------------------------------------------------------------
% Interpolation lineaire de l'efficacite pour determiner la probabilite qu'un evt a d'etre garde
%-----------------------------------------------------------------------------------------------

% M = length(eff_unblend.logE_min);
% te_graph = [eff_unblend.log_tE_min(1) ; eff_unblend.log_tE_min ; eff_unblend.log_tE_max; eff_unblend.log_tE_max(M)];
% eff_graph = [0 ; eff_unblend.efficiency(sort([1:M 1:M])) ; 0];
% figure(1)
% loglog(10.^sort(te_graph), eff_graph)

%-----------------------------------------
%Efficacité en créneau
%---------------------------------
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

%-----------------------------------------
%Efficacité avec interpolation
%---------------------------------
% te_inter = 10.^(eff_field.log_tE_min);
% 
% teffmaxm=max(te_inter);
% teffminm=min(te_inter);
% 
% i1_unblend = find((te<=teffmaxm)&(te>=teffminm));
% i1_blend = find((teblend<=teffmaxm)&(teblend>=teffminm));
% 
% eff_unblend = zeros(1,length(te));	% applique une efficacite nulle aux durees superieures et inferieures
% eff_blend = zeros(1,length(teblend));
% 
% eff_unblend(i1_unblend) = interp1(te_inter,eff_field.efficiency,te(i1_unblend));
% eff_blend(i1_blend) = interp1(te_inter,eff_field.efficiency,teblend(i1_blend));

%--------------------------------------------------------------------------------------------------------------------------
% compare le nombre aleatoire precedent a l'efficacite que l'on vient de calculer afin de decider si l'evt est garde ou non
%--------------------------------------------------------------------------------------------------------------------------

%tirage au sort pour l'efficacité
% ra_unblend = rand(1,length(te))*max(eff_field.efficiency);
% ra_blend = rand(1,length(teblend))*max(eff_field.efficiency);

ra_unblend = rand(1,length(te))*max(eff_field.efficiency);
ra_blend = rand(1,length(teblend))*max(eff_field.efficiency);

% On choisit l'efficacité ici en prenant les bons indices i

i = find(ra_unblend-eff_unblend<=0); 
teobs = te(i);

ib = find(ra_blend-eff_blend<=0); 
teobsblend = teblend(ib); % On récupère les éléments qui sont soumis au blending avec le calcul d'avant

%--------------------
%Calcul gamma avec données éxpérimentale
%-------------------------------
gam_ogle = 0;
gam_obs_ogle = 0;
gam_obs_star_count = 0;
expo = zeros(size(id_field));
% for i = 1:length(id_field)
% % for i = 1:10
%     id_eve = id_field(i);
%     expo(i) = 2741*table5.N21(table5.field == table3.field(id_eve)) / (365.25 * 1e6);
%     gam_obs_star_count = gam_obs_star_count + table3.weight(id_eve)/expo(i);
% %     disp([num2str(gam_obs_star_count), ' ', num2str(expo_i)]);
% 
% end
disp(['gama avec star count bon ', num2str(gam_obs_star_count)])
for i = 1:length(te_inter_min)
    i1_unblend = find(teobs>=te_inter_min(i) & teobs<=te_inter_max(i));
    gam_ogle = gam_ogle + length(i1_unblend)/eff_field.efficiency(i);
    
    i1 = find(teff>=te_inter_min(i) & teff<=te_inter_max(i));
    gam_obs_ogle = gam_obs_ogle + length(i1)/eff_field.efficiency(i);
end

disp(['gamma obs calculé comme ogle',num2str(gam_ogle/length(teobs) / exposure)])
disp(['gamma obs calculé comme ogle avec les donnée ogle efficiency= ',num2str(gam_obs_ogle / exposure)])
disp(['gamma obs calculé comme ogle avec les donnée ogle efficiency= ',num2str(sum(table3.weight(id_field)) / exposure)])
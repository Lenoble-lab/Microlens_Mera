%---------------------------------------------------------------------------------------------------------------
%Microlensing optical depth and event rate in the OGLE-IV Galactic plane fields
%
% P. Mroz, A. Udalski, M.K. Szymanski, I. Soszynski, P. Pietrukowicz,
% S. Kozlowski, J. Skowron, R. Poleski, K. Ulaczyk, M. Gromadzki,
% K. Rybicki, P. Iwanek, and M. Wrona
%---------------------------------------------------------------------------------------------------------------
close all
clear
% Exposure ?


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


%-----------------
%Table 3. Best-fitting parameters of the analyzed microlensing events in low-cadence OGLE fields.
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
table3 = readtable('../OGLEIV/table3.dat',opts);

teff = table3.tE_best;

nbre_bin =  300;
bin_max = 300;
[hist_teff, edges] = histcounts(teff, nbre_bin, 'BinLimits',[0,bin_max], 'BinMethod', 'sturges');

centre = zeros(size(edges)-[0,1]);
for j =1:length(centre);
centre(j)=(edges(j)+edges(j+1))/2;
end

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
long = 1;

%exp ogle
i0 = find(abs(table6.glon - long)<0.5 & table6.glat<0);

%Calcul modèle
tau_load = load('graph_iso_model.mat');

[L, B] = meshgrid(tau_load.L, tau_load.B);

i1 = find(tau_load.L==long);
i2 = find(tau_load.B<0);

figure(1)
hold on
errorbar(table6.glat(i0), table7.tau(i0), table7.tau_err(i0), 'o')
plot(tau_load.B(i2), tau_load.tau_table(i1,i2)*1e6)
legend('Mesure d''OGLE IV', 'Modèle')
xlabel('b (deg)')
ylabel('\tau \times 10^{-6}')
%%
%---------------
%Efficacité
%-----------
VarNames = {'log_tE_min', 'log_tE_max', 'efficiency'};
VarType = {'double', 'double', 'double'};

opts = delimitedTextImportOptions('VariableNames',VarNames,'VariableTypes',VarType,...
        'Delimiter',' ', 'DataLines', 4, 'WhiteSpace', ' ', 'ConsecutiveDelimitersRule', 'join');
    
% eff = readtable(strcat('../OGLEIV/eff/',extractBetween(table3.field(1), 1, 6),'.eff'),opts);
eff = readtable(strcat('../OGLEIV/eff/','BLG500','.eff'),opts);

%-----------------------------------------------------------------------------------------------
% Interpolation lineaire de l'efficacite pour determiner la probabilite qu'un evt a d'etre garde
%-----------------------------------------------------------------------------------------------

figure(1)
plot((10.^eff.log_tE_min + 10.^eff.log_tE_max)/2, eff.efficiency)


[tinterpmacho,indice] = sort(tmacho2005); 
effinterpmacho = effmacho2005t(indice) ;

%Il y a des doublons dans les données, il faut les supprimer pour que l'interpolation se passe correctement
indices = [1:length(tinterpmacho)-1];
il = find(tinterpmacho(indices)~=tinterpmacho(indices+1));
tinterpmacho = [tinterpmacho(il),tinterpmacho(length(tinterpmacho))];
effinterpmacho = [effinterpmacho(il),effinterpmacho(length(tinterpmacho))];


teffmaxm=max(tinterpmacho);
teffminm=min(tinterpmacho);

i1 = find((te<=teffmaxm)&(te>=teffminm));
effsimmacho = zeros(1,length(te));	% applique une efficacite nulle aux durees superieures et inferieures
effsimmachoblend = zeros(1,length(te));
effsimmacho(i1) = interp1(tinterpmacho,effinterpmacho,te(i1));
effsimmachoblend(i1) = interp1(tinterpmacho,effinterpmacho,teblend(i1));

%--------------------------------------------------------------------------------------------------------------------------
% compare le nombre aleatoire precedent a l'efficacite que l'on vient de calculer afin de decider si l'evt est garde ou non
%--------------------------------------------------------------------------------------------------------------------------

%tirage au sort pour l'efficacité
ramacho = rand(1,length(te))*max(effinterpmacho);

% On choisit l'efficacité ici en prenant les bons indices i

i = find(ramacho-effsimmacho<=0); 

teobs = te(i);

ib = find(ramacho-effsimmachoblend<=0); 

teobsblend = teblend(ib); % On récupère les éléments qui sont soumis au blending avec le calcul d'avant


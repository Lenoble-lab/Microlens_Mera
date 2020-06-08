%---------------------------------------------------------------------------------------------------------------
%OGLE-III MICROLENSING EVENTS AND THE STRUCTURE OF THE GALACTIC BULGE ∗
% Łukasz Wyrzykowski 1,2,5 , Alicja E. Rynkiewicz 1 , Jan Skowron 1 , Szymon Kozłlowski 1 , Andrzej Udalski 1 ,
% Michał
% l K. Szymański 1 , Marcin Kubiak 1 , Igor Soszyński 1 , Grzegorz Pietrzyński 1,3 , Radosłlaw Poleski 1,4 ,
% Paweł
% l Pietrukowicz 1 , and Michałl Pawlak
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
%Table 2 : Microlensing parameter for 3500 standard OGLE III
% Microlensing events of class A
%------------------------------------------------
delimiter = ' ';
VarNames_table2 = {'id', 'ra', 'dec', 'field', 'N_stars', 't0', 'e_t0', 'E_t0', 'te', 'e_te', 'E_te', ...
'u0', 'e_u0', 'E_u0', 'fs', 'e_fs', 'E_fs', 'I0', 'e_I0', 'E_I0', 'Chi', 'Ndof', 'EWS_id', 'glon', 'glat'};
VarTypes_table2 = {'double', 'string', 'string', 'double', 'double', 'double', 'double', ...
'double', 'double', 'double', 'double','double', 'double', 'double', 'double', 'string'}; 

opts = delimitedTextImportOptions('VariableNames',VarNames_table2,'VariableTypes',VarTypes_table2,...
                                'Delimiter',delimiter, 'DataLines', 8, ...
                       'WhiteSpace', ' ', 'ConsecutiveDelimitersRule', 'join');
table2 = readtable('../OGLEIII/table2.table',opts);



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


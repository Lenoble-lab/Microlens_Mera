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
'u0', 'e_u0', 'E_u0', 'fs', 'e_fs', 'E_fs', 'I0', 'e_I0', 'E_I0', 'Chi', 'Ndof','V' ,'e_V', 'EWS_id', 'glon', 'glat'};
VarTypes_table2 = {'double', 'string', 'string', 'string', 'double', 'double', 'double', ...
'double', 'double', 'double', 'double','double', 'double', 'double', 'double', ...
'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'string', 'double', 'double'}; 

opts = delimitedTextImportOptions('VariableNames',VarNames_table2,'VariableTypes',VarTypes_table2,...
                                'Delimiter',delimiter, 'DataLines', 8, ...
                       'WhiteSpace', ' ', 'ConsecutiveDelimitersRule', 'join');
table2 = readtable('../OGLEIII/table2.txt',opts);

%remet glon dans [-180, 180]:
re = rem(table2.glon, 180);
i0 = find(table2.glon~=re);
table2.glon(i0) = re(i0) - ones(size(i0)).*180;

%--------------------------
%Efficacité
%---------------------------

VarNames_eff = {'id', 't_e', 'N_events', 'undef', 'efficiency'};
VarTypes_eff = {'double', 'double', 'double', 'double', 'double'};

opts = delimitedTextImportOptions('VariableNames',VarNames_eff,'VariableTypes',VarTypes_eff,...
                                'Delimiter',delimiter,...
                       'WhiteSpace', ' ', 'ConsecutiveDelimitersRule', 'join');
eff = readtable('../OGLEIII/efficiency.dat',opts);

M = 35;
figure
loglog(eff.t_e(sort([1:M 1:M])), [0 ; eff.efficiency(sort([1:M-1 1:M-1])) ; 0], '-')
hold on
loglog(eff.t_e, eff.efficiency)

%---------------
%Tracé fig. 1, pour vérifier la transformation des angles en coordonnées
%galactique
%-------------
% figure(2)
% plot(table2.glon, table2.glat, 'x')

%----------------------
%tracé figure 7
%----------------------

%Choix des évènements avec une erreur relative raisonable
rel_err = (abs(table2.e_te) + table2.E_te)./(2.*table2.te);
i0 = find(rel_err<0.5);

figure(2)
semilogxhistnormalise(table2.te(i0), 25)




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
%-----------------------------------------------------------------------------------------------
% Interpolation lineaire de l'efficacite pour determiner la probabilite qu'un evt a d'etre garde
%-----------------------------------------------------------------------------------------------

figure(1)
plot((10.^eff.log_tE_min + 10.^eff.log_tE_max)/2, eff.efficiency)


[tinterpmacho,indice] = sort(tmacho2005); 
effinterpmacho = effmach
o2005t(indice) ;

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

function semilogxhistnormalise(val,M)
% semilogxhist - generate histogram with M bars and log-scale x axis
vmin=min(val); vmax=max(val);
edges=vmin*(vmax/vmin).^([0:M]/M);
count=histcounts(val,edges, 'Normalization', 'probability'); 
if size(count,2)==1, count=count'; end 
x=edges(sort([1:M 1:M])); 
y=[0 count(sort([1:M-1 1:M-1])) 0];

% outline only: semilogx(x, y, '-');
plot(x, y, '-'); 
% fill(x, y, 'b'); 
set(gca,'XScale','log');
end
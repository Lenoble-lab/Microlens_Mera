%---------------------------------------------------------------------------------------------------------------
%OGLE-III MICROLENSING EVENTS AND THE STRUCTURE OF THE GALACTIC BULGE ∗
% Łukasz Wyrzykowski 1,2,5 , Alicja E. Rynkiewicz 1 , Jan Skowron 1 , Szymon Kozłlowski 1 , Andrzej Udalski 1 ,
% Michał
% l K. Szymański 1 , Marcin Kubiak 1 , Igor Soszyński 1 , Grzegorz Pietrzyński 1,3 , Radosłlaw Poleski 1,4 ,
% Paweł
% l Pietrukowicz 1 , and Michałl Pawlak
%---------------------------------------------------------------------------------------------------------------
% close all
% clear
% Exposure ?
exposure = 400;

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
i_err = find(table2.glon~=re);
table2.glon(i_err) = re(i_err) - ones(size(i_err)).*180;

teff = table2.te;
%--------------------------
%Efficacité
%---------------------------

VarNames_eff = {'id', 't_e', 'N_events', 'undef', 'efficiency'};
VarTypes_eff = {'double', 'double', 'double', 'double', 'double'};

opts = delimitedTextImportOptions('VariableNames',VarNames_eff,'VariableTypes',VarTypes_eff,...
                                'Delimiter',delimiter,...
                       'WhiteSpace', ' ', 'ConsecutiveDelimitersRule', 'join');
eff_field = readtable('../OGLEIII/efficiency.dat',opts);

M = 35;
figure(1)
loglog(eff_field.t_e(sort([1:M 1:M])), [0 ; eff_field.efficiency(sort([1:M-1 1:M-1])) ; 0], '-')
hold on
loglog(eff_field.t_e, eff_field.efficiency)

eff = eff_field.efficiency;

%----------------------
%tracé figure 7 (histogramme des évènements)
%----------------------

% %Choix des évènements avec une erreur relative raisonable
rel_err = max(abs(table2.e_te), table2.E_te)./(table2.te);
i_err = find(rel_err<0.5);
i_BW_dir = find(abs(table2.glon-1)<2 & abs(table2.glat+2)<2);
i_BW = intersect(i_BW_dir, i_err);


% 
% figure(2)
% semilogxhistnormalise(table2.te(i0), 25)
% [hist_teff, edges] = histcounts(table2.te, nbre_bin, 'BinLimits',[0,bin_max], 'BinMethod', 'sturges');
% 
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


%--------------------------
%Figure 13 (histogramme corrigé de l'efficacité)
%---------------------
bin_max = 300;
nbre_bin = 20;
v = table2.te;
weight = zeros(size(table2.te)); %efficacité

%table d'efficacité pour chaque événements
for i = 1:length(weight)
    i0 = find(v(i)<eff_field.t_e,1);
    weight(i) = 1/eff_field.efficiency(i0);
end

[histw, vinterval] = histwc(table2.te(i_err), weight(i_err), 50, 1, 400);
[histw_BW, vinterval_BW] = histwc(table2.te(i_BW), weight(i_BW), 50, 1, 400);

M = length(vinterval);
M_BW = length(vinterval_BW);
figure(3)
loglog(vinterval(sort([1:M 1:M 1:1] )), [0; 0 ; histw(sort([1:M-1 1:M-1])) ; 0])
hold on
loglog(vinterval_BW(sort([1:M 1:M 1:1] )), [0; 0 ; histw_BW(sort([1:M-1 1:M-1])) ; 0])

%-----------------------------------------
%Efficacité avec interpolation
%---------------------------------

teffmaxm=max(eff_field.t_e);
teffminm=min(eff_field.t_e);
i1_unblend = find((te<=teffmaxm)&(te>=teffminm));
i1_blend = find((teblend<=teffmaxm)&(teblend>=teffminm));

eff_unblend = zeros(1,length(te));	% applique une efficacite nulle aux durees superieures et inferieures
eff_blend = zeros(1,length(teblend));
% 
eff_unblend(i1_unblend) = interp1(eff_field.t_e,eff_field.efficiency,te(i1_unblend));
eff_blend(i1_blend) = interp1(eff_field.t_e,eff_field.efficiency,teblend(i1_blend));

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

function [histw, vinterval] = histwc(vv, ww, nbins, minn, maxx)
  minV  = minn;
  maxV  = maxx;
  delta = (maxV-minV)/nbins;
  vinterval = 10.^linspace(log(minV), log(maxV), nbins);
  histw = zeros(nbins, 1);
  for i=1:length(vv)
    ind = find(vinterval < vv(i), 1, 'last' );
    if ~isempty(ind)
      histw(ind) = histw(ind) + ww(i);
    end
  end
end

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
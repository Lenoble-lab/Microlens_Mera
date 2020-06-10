%---------------------------------------------------------------------------------------------------------------
%POSSIBLE SOLUTION OF THE LONG-STANDING DISCREPANCY IN THE MICROLENSING OPTICAL DEPTHTOWARD 
%THE GALACTIC BULGE BY CORRECTING THE STELLAR NUMBER COUNTT. 
%Sumi and M. T. Penny
%AUssi tiré de Sumi et al 2013 (MOA II)
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
%Table 4. Average microlensing optical depth and event rates at the position 
%of each subfield for the all-source sample
%------------------------------------------------

delimiter = ' ';
VarNames_table7_MOA = {'blank', 'Field', 'glon', 'glat',  'Nsub',  'Nstar',  'Nev', 'tau',  'etau', 'e_tau', 'E_tau', 'Gamma', 'eGamma', ...
'e_Gamma', 'E_Gamma', 'Gammad', 'eGammad', 'e_Gammad', 'E_Gammad'};
VarTypes_table4_MOA = {'string', 'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double',...
 'double', 'double', 'double', 'double', 'double', 'double', 'double'}; 

opts = delimitedTextImportOptions('VariableNames',VarNames_table7_MOA,'VariableTypes',VarTypes_table4_MOA,...
                                'Delimiter',delimiter, 'DataLines', 22, ...
                       'WhiteSpace', '  ', 'ConsecutiveDelimitersRule', 'join');
table7_MOA = readtable('../MOA_II/Table7.dat',opts);

%---------------------------------------------------
%Table 4 : Microlensing events used in the optical depth and event rate measurements.
%-----------------------------------------------------
VarNames_table4_MOA = {'blank',  'ID', 'RA', 'Dec', 'Ndata', 't0', 'tE', 'e_tE', 'u0', 'e_u0', 'Is', 'chi2dof'};
VarTypes_table4_MOA = {'string', 'string', 'string', 'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'};

opts = delimitedTextImportOptions('VariableNames',VarNames_table4_MOA,'VariableTypes',VarTypes_table4_MOA,...
                                'Delimiter',delimiter, 'DataLines', 22, ...
                       'WhiteSpace', '  ', 'ConsecutiveDelimitersRule', 'join');
                   
table4_MOA = readtable('../MOA_II/Table4.dat',opts);





%----------------
%Tracé profondeur optique en fonction de la lattitude galactique au centre
%---------------------
% glong = 0;
% glat = (0:.1:1) .*-6 -1;
% 
% tau_mean = zeros(size(glat)-1);
% lat_mean = zeros(size(glat)-1);
% e_tau_mean = zeros(size(glat)-1);
% E_tau_mean = zeros(size(glat)-1);
% 
% for i = 2:length(glat)
% i0 = find(abs(table7_MOA.glon - glong)<5 & glat(i-1)>table7_MOA.glat & table7_MOA.glat>glat(i));
% i = i-1;
% tau_mean(i) = mean(table7_MOA.tau(i0));
% lat_mean(i) = mean(table7_MOA.glat(i0));
% e_tau_mean(i) = mean(table7_MOA.e_tau(i0));
% E_tau_mean(i) = mean(table7_MOA.E_tau(i0));
% 
% disp([num2str(lat_mean(i)), '    ', num2str(mean(table7_MOA.tau(i0)))])
% disp(['nbre subfield ', num2str(length(i0))])
% disp(i)
% end

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


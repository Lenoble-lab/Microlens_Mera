%---------------------------------------------------------------------------------------------------------------
%POSSIBLE SOLUTION OF THE LONG-STANDING DISCREPANCY IN THE MICROLENSING OPTICAL DEPTHTOWARD 
%THE GALACTIC BULGE BY CORRECTING THE STELLAR NUMBER COUNTT. 
%Sumi and M. T. Penny
%AUssi tiré de Sumi et al 2013 (MOA II)
%---------------------------------------------------------------------------------------------------------------

% close all
% clear
% Exposure 
exposure = 569/365.25 * 110.3;

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

%-----------------------------------------
%Pour l'instant, pas trouvé l'efficacité de l'exp alors prend eff de 1
%--------------------------------------------

teobs = te;

teobsblend = teblend; % On récupère les éléments qui sont soumis au blending avec le calcul d'avant


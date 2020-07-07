% Programme de micro-lentilles gravitationnelles (bulbe et disque)

clearvars -except field field_li

%---------------------------------------------------------------------------------------------------------------
%Microlensing optical depth and event rate in the OGLE-IV Galactic plane fields
%
% P. Mroz, A. Udalski, M.K. Szymanski, I. Soszynski, P. Pietrukowicz,
% S. Kozlowski, J. Skowron, R. Poleski, K. Ulaczyk, M. Gromadzki,
% K. Rybicki, P. Iwanek, and M. Wrona
%---------------------------------------------------------------------------------------------------------------

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
eff_field = readtable('../OGLEIV/table5_central.txt',opts);

field_list = ["BLG500", "BLG501", "BLG504", "BLG505", "BLG506", "BLG511", "BLG512", "BLG534", "BLG611"];

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

%-----------------------------------------------
%Raw te
%-----------------------------------------------

VarNames_table_event = {'Name', 't_0', 't_E', 'u_0'};
VarTypes_table_event = {'string', 'double', 'double', 'double'}; 

opts = delimitedTextImportOptions('VariableNames',VarNames_table_event,'VariableTypes',VarTypes_table_event,...
                                'Delimiter',delimiter, 'DataLines', 8, ...
                       'WhiteSpace', ' ', 'ConsecutiveDelimitersRule', 'join');
table_event = readtable('../OGLEIV/OGLE-IV-events-FFP.txt',opts);


%-------------
% Choix du champ à analyser
%------------

field = "BLG505";

%-------------
%donnée expérimentales
%----------------

i_field = find(extractBetween(table_event.Name, 1,6) == field);
teff = table_event.t_E(i_field);

%------------------
%fichiers resultats
%------------------

fichevents='evenements.sim1.txt';

global vsr vsp vst vlp vlt vlr
%----------------------------------
% Constantes physiques (unites SI)
%----------------------------------
G=6.672e-11;	pc=3.08567802e16;
kpc=pc*1e3;  	Msol=1.989e30;
c=299792458;	GMsol=1.32712497e20;

%-------------------------------------------------
% Param�tres vlimit (vitesse perpendiculaire maxi)
%-------------------------------------------------

global vlimit

vlimit = 1000e3;

%----------------------
% Nombre de simulations
%----------------------

n = 20000;
% n = 5000;
nbsimul=500; %a augmenter pour meilleure stat
nbMAX=500;

%----------------------------------------------------------------------
% Param�tres de la fonction de distribution de la distance de la source
%----------------------------------------------------------------------

global dsup dinf 

dsup = 15000.;
dinf = 800.;
%distance en parsec

%-------
% soleil 
%-------
global Ro elev
Ro = 8000;	   % distance Sun-GC en pc
elev = 26 ;      % elevation au dessus du plan du disque
Lkpc=Ro/1000;  % distance Sun-GC en kpc


%----------------------------
%rayon de corotation (en pc)
%----------------------------

global Rcoro
Rcoro = 3500;


%---------------------------------------       
% param�tres du calcul de microlentilles
%---------------------------------------
global l b

% on suit la direction du champ considéré

l = table6.glon(table6.field == field) *pi/180;    % direction d'observation en radian
b = table6.glat(table6.field == field) *pi/180;

uT = 1;		   % Seuil de d�tection en param�tre d'impact
AT = 3/sqrt(5);    % Seuil de d�tection en amplification


%-----------------------------------
% Param�tres de la fonction de masse
%-----------------------------------
global minfbu msupbu minfdm msupdm minfde msupde minfh msuph 
global minf msup

minf=0.01;
msup=100;

    %--------------------
    % cas du disque mince
    %--------------------

minfdm = minf;
msupdm = msup;

    %--------------------
    % cas du disque epais
    %--------------------

minfde = minf;
msupde = msup;

    %-------------
    % cas du bulbe
    %-------------

minfbu = minf;
msupbu = msup;

    %------------
    % cas du halo
    %------------

minfh = minf;
msuph = msup;


msuptot=max([msupdm msupde msupbu msuph]);
minftot=max([minfdm minfde minfbu minfh]);

%-------------------------------------------------
% Calcul pr�alable du cosinus pour aller plus vite
%-------------------------------------------------
global sinb cosb  cosbl sinl cosl
sinb = abs(sin(b));		cosb = cos(b);		cosl = cos(l);
cosbl=cos(b)*cos(l);		sinl = sin(l);


%----------------------------------------
% Param�tres de la fonction de luminosité pour le blending (Alibert et al.)
%----------------------------------------

global Vinf Vsup normfl

Vinf = 22; % magnitudes limites des étoiles observées.
Vsup = 16;

Linf=lum(Vinf);
Lsup=lum(Vsup);

% Norme de la fonction de luminosité :

normfl = integral(@fl, Linf, Lsup);


%-------------------------------
%calcul de la profondeur optique
%-------------------------------

normnu=integral(@nsource,dinf,dsup);
Ctau = 4*pi*GMsol*uT*uT/c/c/pc;
tau = Ctau * integral2(@dtau,0,1,dinf,dsup, 'Method', 'iterated') /normnu;
taucx=tau;
tau=real(taucx);
disp(['tau = ' num2str(taucx)]);

%-----------------------------------------------------------------------------
%calcul de la profondeur optique dans le cas d'une source a distance constante
%-----------------------------------------------------------------------------

%Ctau = 4*pi*GMsol*uT*uT/c/c/pc;
%tau  = Ctau*quad('dtaud',0,1,1e-2);
%disp(['tau = ' num2str(tau)]);

%----------------
%liste pour observer ensuite la pdmf
%----------------------
m_tot_imf = [];
m_tot_pdmf = [];
frac_M_tot = [];
frac_N_tot = [];

%-------------------------------------------------------
%calcul de la masse moyenne et de la normalisation de fm
%-------------------------------------------------------

    %--------------------
    % cas du disque mince
    %--------------------

global normfmdm mmeandm

normfmdm=integral(@fmdm,minfdm,msupdm);
disp(['integrale de la fonction de masse du disque mince = ' num2str(normfmdm)]);

mmeandm=integral(@mPmdm,minfdm,msupdm);

    %--------------------
    % cas du disque epais
    %--------------------

global normfmde mmeande

normfmde=integral(@fmde,minfde,msupde);
disp(['integrale de la fonction de masse du disque epais = ' num2str(normfmde)]);

mmeande=integral(@mPmde,minfde,msupde);


    %-------------
    % cas du bulbe
    %-------------

global normfmbu mmeanbu

normfmbu=integral(@fmbu,minfbu,msupbu);
disp(['integrale de la fonction de masse du bulbe = ' num2str(normfmbu)]);

mmeanbu=integral(@mPmbu,minfbu,msupbu);

    %------------
    % cas du halo
    %------------

    global normfmh mmeanh

normfmh=integral(@fmh,minfh,msuph);
disp(['integrale de la fonction de masse du halo = ' num2str(normfmh)]);

mmeanh=integral(@mPmh,minfh,msuph);

%-------------------------------------------------------------
% Preparation de la table pour le tirage aleatoire de la masse
%-------------------------------------------------------------

    %--------------------
    % cas du disque mince
    %--------------------

mmdm = (0:1e-5:1).*(msupdm-minfdm)+minfdm;
fmmdm = probadm(mmdm);
ifmmdm = cumsum(fmmdm);	% primitive de la fonction de masse fm
ifmmdm = ifmmdm-ifmmdm(1);
ifmmdm = ifmmdm./ifmmdm(end);	% on fait en sorte que la primitive varie de 0 a 1

    %--------------------
    % cas du disque epais
    %--------------------

mmde = (0:1e-5:1).*(msupde-minfde)+minfde;
fmmde = probade(mmde);
ifmmde = cumsum(fmmde);	% primitive de la fonction de masse fm
ifmmde = ifmmde-ifmmde(1);
ifmmde = ifmmde./ifmmde(end);	% on fait en sorte que la primitive varie de 0 a 1

    %-------------
    % cas du bulbe
    %-------------
 
mmbu = (0:1e-5:1).*(msupbu-minfbu)+minfbu;
fmmbu = probabu(mmbu);
ifmmbu = cumsum(fmmbu);	% primitive bu la fonction bu masse fm
ifmmbu = ifmmbu-ifmmbu(1);
ifmmbu = ifmmbu./ifmmbu(end);	% on fait en sorte que la primitive varie de 0 a 1

    %------------
    % cas du halo
    %------------

mmh = (0:1e-5:1).*(msuph-minfh)+minfh;
fmmh = probah(mmh);
ifmmh = cumsum(fmmh);	% primitive h la fonction h masse fm
ifmmh = ifmmh-ifmmh(1);
ifmmh = ifmmh./ifmmh(end);	% on fait en sorte que la primitive varie de 0 a 1

%-----------------------------------------------------------------------------
% Preparation de la table pour le tirage aleatoire de la distance de la source
%-----------------------------------------------------------------------------

% dd = (0:1e-5:1).*(dsup-dinf)+dinf;
% fds = nsource (dd);
% ifds = real(cumsum(fds));	% primitive 
% ifds = ifds-ifds(1);
% ifds = ifds./ifds(end);	% on fait en sorte que la primitive varie de 0 a 1

dd = (0:1e-5:1).*(dsup-dinf)+dinf;
disp(dd(500));
fds = nsource (dd);
ifds = cumsum(fds);	% primitive 
ifds = ifds-ifds(1);
ifds = real(ifds./ifds(end));	% on fait en sorte que la primitive varie de 0 a 1


%-------------------------------------------
% calcul du maximum de Gamma par monte carlo
%-------------------------------------------

for compteur = 1:nbMAX,
    
    if ( rem(compteur,100) == 0 )
    
disp(['compteur = ' num2str(compteur)]);

    end

% randn('seed',sum(100000*clock)), rand('seed',sum(100000*clock))

x=rand(1,n);
ds=rand(1,n)*(dsup-dinf)+dinf;
m=rand(1,n)*(msuptot-minftot)+minftot;
v=rand(1,n)*vlimit;

gmax(compteur)=max(dgam(x,ds,v,m));

clear x ds m v

end

Gammax=max(gmax);
disp(['Gammax = ' num2str(Gammax)]);


%---------------------
%---------------------
% debut du Monte-Carlo
%---------------------
%---------------------


%----------------------------------
% ouverture des fichiers de donnees
%----------------------------------

fid = fopen(fichevents,'w');
fid1 = fopen('evenements.txt','w');
fids = fopen('strange.txt','w');

dgamma=[];
dgammaccepte=[];        % initialisation pour les traitement ulterieur des donnees
tecorrespondant=[];

for compteur = 1:nbsimul,

  if ( rem(compteur,100) == 0 )  
    
disp(['compteur = ' num2str(compteur)]);

  end


%---------------------------------------------------
% Initialisation aleatoire des generateurs aleatoire
%---------------------------------------------------

% randn('seed',sum(100000*clock)), rand('seed',sum(100000*clock))

%----------------------------------
%tirage de la distance de la source
%----------------------------------

ra = rand(1,n);
[ifds, index] = unique(ifds); 

ds = interp1(ifds,dd(index),ra);

%------------
% Tirage de x
%------------

x=rand(1,n);

%------------------------------------------
%a quelle population appartient la source ?
%------------------------------------------

ra=rand(1,n);
[R,z,th]= toGC(ds);
rhotot = rhodm(R,z,th) + rhode(R,z,th) + rhobulbe(R,z,th) + rhohalo(R,z,th);
idms=find(ra<=(rhodm(R,z,th)./rhotot));
ides=find(rhodm(R,z,th)./rhotot<ra & ra <= (rhodm(R,z,th)+rhode(R,z,th))./rhotot);
ibus=find((rhodm(R,z,th)+rhode(R,z,th))./rhotot<ra & ra <= (rhodm(R,z,th)+rhode(R,z,th)+rhobulbe(R,z,th))./rhotot);
ihs=find(ra >= (rhodm(R,z,th)+rhode(R,z,th)+rhobulbe(R,z,th))./rhotot);

%--------------------------------------------
%a quelle population appartient la lentille ?
%--------------------------------------------

ra=rand(1,n);
[R,z,th]= toGC(x.*ds);
rhotot = rhodm(R,z,th) + rhode(R,z,th) + rhobulbe(R,z,th) + rhohalo(R,z,th);
idml=find(ra<=(rhodm(R,z,th)./rhotot));
idel=find(rhodm(R,z,th)./rhotot<ra & ra <= (rhodm(R,z,th)+rhode(R,z,th))./rhotot);
ibul=find((rhodm(R,z,th)+rhode(R,z,th))./rhotot<ra & ra <= (rhodm(R,z,th)+rhode(R,z,th)+rhobulbe(R,z,th))./rhotot);
ihl=find(ra >= (rhodm(R,z,th)+rhode(R,z,th)+rhobulbe(R,z,th))./rhotot);


%-----------------------------------------------------------------
% test : methode de Peale. Lentilles et sources ont les proprietes
% du disque si R>Rcoro, et du bulbe sinon
%-----------------------------------------------------------------

% [R,z,th]= toGC(x.*ds);
% ibul=find(R<=Rcoro);
% idml=find(R>Rcoro);
% idel=find(R<0.);
% ihl=find(R<0.);

% [R,z,th]= toGC(ds);
% ibus=find(R<=Rcoro);
% idms=find(R>Rcoro);
% ides=find(R<0.);
% ihs=find(R<0.);


%---------------------------------
%tirage de la masse de la lentille
%---------------------------------

ra=rand(1,n);
[ifmmdm, index] = unique(ifmmdm); 
m(idml) = interp1(ifmmdm,mmdm(index),ra(idml));
[ifmmde, index] = unique(ifmmde); 
m(idel) = interp1(ifmmde,mmde(index),ra(idel));
[ifmmbu, index] = unique(ifmmbu); 
m(ibul) = interp1(ifmmbu,mmbu(index),ra(ibul));
[ifmmh, index] = unique(ifmmh); 
m(ihl) = interp1(ifmmh,mmh(index),ra(ihl));

%-----------------
%PDMF
%------------------

m_tot_imf = [m_tot_imf, m];
% [m, frac_N, frac_M] = PDMF_gould_1(m);
[m, frac_N, frac_M] = PDMF_Maraston_1(m);

m_tot_pdmf = [m_tot_pdmf, m];
frac_N_tot = [frac_N_tot ; frac_N];
frac_M_tot = [frac_M_tot ; frac_M];

% m = zeros(size(x)) + 0.6;
%----------------------------------------------
%calcul des sigmas et v rotation de la lentille
%----------------------------------------------
[R,z,th]= toGC(x.*ds);  % Attention  car il faut recalculer R qui a ete change pour le calcul toGC pour la source


%direction u_r et u_theta

%Disque mince
sigrl(idml)=sigrdm(R(idml),z(idml),th(idml));
sigtl(idml)=sigtdm(R(idml),z(idml),th(idml));

%Disque épais
sigrl(idel)=sigrde(R(idel),z(idel),th(idel));
sigtl(idel)=sigtde(R(idel),z(idel),th(idel));

%Bulbe
sigrl(ibul)=sigrb(R(ibul),z(ibul),th(ibul));
sigtl(ibul)=sigtb(R(ibul),z(ibul),th(ibul));

%Halo
% sigrl(ihl)=sigrh(R(ihl),z(ihl),th(ihl));
% sigtl(ihl)=sigth(R(ihl),z(ihl),th(ihl));


%Coordonnée cylindrique u_z
sigzl(idml)=sigzdm(R(idml),z(idml),th(idml));
sigzl(ibul)=sigzb(R(ibul),z(ibul),th(ibul));
sigzl(idel)=sigzde(R(idel),z(idel),th(idel));
% sigzl(ihl)=sigzh(R(ihl),z(ihl),th(ihl));

%Coordonnée sphérique u_phi
% sigpl(idml)=sigpdm(R(idml),z(idml),th(idml));
% sigpl(idel)=sigpde(R(idel),z(idel),th(idel));
% sigpl(ibul)=sigpb(R(ibul),z(ibul),th(ibul));
% sigpl(ihl)=sigph(R(ihl),z(ihl),th(ihl));


vrotl(idml)=vrotdm(R(idml),z(idml),th(idml));
vrotl(idel)=vrotde(R(idel),z(idel),th(idel));
vrotl(ibul)=vrotb(R(ibul),z(ibul),th(ibul));
vrotl(ihl)=vroth(R(ihl),z(ihl),th(ihl));

%--------------------------------------------
%calcul des sigmas et v rotation de la source
%--------------------------------------------

[R,z,th]= toGC(ds);  % Attention  car il faut recalculer R

%disque mince
sigrs(idms)=sigrdm(R(idms),z(idms),th(idms));
sigts(idms)=sigtdm(R(idms),z(idms),th(idms));

%disque épais
sigrs(ides)=sigrde(R(ides),z(ides),th(ides));
sigts(ides)=sigtde(R(ides),z(ides),th(ides));

%Bulbe
sigrs(ibus)=sigrb(R(ibus),z(ibus),th(ibus));
sigts(ibus)=sigtb(R(ibus),z(ibus),th(ibus));

%Halo
% sigrs(ihs)=sigrh(R(ihs),z(ihs),th(ihs));
% sigts(ihs)=sigth(R(ihs),z(ihs),th(ihs));

%Coordonnée cylindrique u_z
sigzs(idms)=sigzdm(R(idms),z(idms),th(idms));
sigzs(ides)=sigzde(R(ides),z(ides),th(ides));
sigzs(ibus)=sigzb(R(ibus),z(ibus),th(ibus));
% sigzs(ihs)=sigzh(R(ihs),z(ihs),th(ihs));


%Coordonnée sphérique u_phi
% sigps(idms)=sigpdm(R(idms),z(idms),th(idms));
% sigps(ides)=sigpde(R(ides),z(ides),th(ides));
% sigps(ibus)=sigpb(R(ibus),z(ibus),th(ibus));
% sigps(ihs)=sigph(R(ihs),z(ihs),th(ihs));



vrots(idms)=vrotdm(R(idms),z(idms),th(idms));
vrots(ides)=vrotde(R(ides),z(ides),th(ides));
vrots(ibus)=vrotb(R(ibus),z(ibus),th(ibus));
vrots(ihs)=vroth(R(ihs),z(ihs),th(ihs));

%cas particulier pour comparer avec Peale, la vitesse de rotation de la source est nulle

%vrots=zeros(size(x));

%---------------------------------------------------------------------------------
% cas particulier ou la source est a une distance fixe, et immobile ( cas du LMC )
%---------------------------------------------------------------------------------

%Lsource=8000; 				 % distance de la source, en parsecs
%ds=Lsource*ones(size(m));
%sigrs=zeros(size(m));
%sigps=zeros(size(m));
%sigts=zeros(size(m));
%vrots=zeros(size(m));

%---------------------------------------------------
% Tirage pour les vitesses, (l)entille puis (s)ource
%---------------------------------------------------

glr = rand(1,n);	gsr = rand(1,n);
glt = rand(1,n);	gst = rand(1,n);
glp = rand(1,n);	gsp = rand(1,n);

%cooronnées phériques
% v = vperp(x,glr,glt,glp,sigrl,sigtl,sigpl,vrotl,ds,gsr,gst,gsp,sigrs,sigts,sigps,vrots);


%coordonées cylindriques
v = vperp_cyl(x,glr,glp,glt,sigrl,sigtl,sigzl,vrotl,ds,gsr,gsp,gst,sigrs,sigts,sigzs,vrots);

vlim=vlimit*ones(size(v));

vi = find(v-vlim>0);
if (length(vi)>=1)
  disp(['vitesse superieure a vlimite : ' num2str(max(v(vi)))]);
end;  

%-------------------------------------------------------------------------
% Tirage du nombre qui permet de dire si un point simule est retenu ou non
%-------------------------------------------------------------------------

ra = rand(size(x));
g  = dgam(x,ds,v,m)./Gammax;
i = find(g-ra>=0);

j = find(g>1);

dgamma = [dgamma,g];

dgammaccepte=[dgammaccepte,g(i)];

%------------------------
%selection des evenements
%------------------------

if (length(i)>=1)
  te = (2/c*sqrt(GMsol*pc)/86400).*sqrt(ds(i).*m(i).*x(i).*(1-x(i)))./v(i);
  
  tecorrespondant=[tecorrespondant,te];
  
  evnts = [x(i);ds(i);v(i);m(i);te];
  fprintf(fid,'%12.8f  %12.8f  %12.8f  %12.8f  %12.8f\n',evnts);
  fprintf(fid1,'%12.8f  %12.8f  %12.8f  %12.8f  %12.8f\n',evnts);
end;


if (length(j)>=1)
  strange = [x(j);ds(j);v(j);m(j)];
  fprintf(fids,'%12.8f  %12.8f  %12.8f  %12.8f\n',strange);
end;


clear x glr glt glp sigrl sigtl sigpl vrotl ds gsr gst gsp sigrs sigts sigps vrots 
clear v m ibul idel idml ihl ibus ides idms ihs ra R z th rhotot evnts te strange

end



fclose(fid);
fclose(fids);


%---------------------
%% Analyse des résultats
%---------------------

close all;
fid=fopen('simul_para.txt','w');

fprintf(fid,'%15.11f  \n',tau)

fprintf(fid,'%15.11f  \n',n);
fprintf(fid,'%15.11f  \n',nbsimul);
fprintf(fid,'%15.11f  \n',tau);
fprintf(fid,'%15.11f  \n',Gammax);
fprintf(fid,'%15.11f  \n',uT);
fprintf(fid,'%15.11f  \n',AT);

fclose(fid);

%-------------------
%-------------------
% fin du Monte-Carlo
%-------------------
%-------------------
%%

%----------------------------------------
% recuperation des evenements selectionnes
%----------------------------------------

load evenements.txt
x=evenements(:,1);
ds=evenements(:,2);
v=evenements(:,3);
m=evenements(:,4);
te=evenements(:,5);

x=x';
ds=ds';
v=v';
m=m';
te=te';

disp(' ')
close all

%------------------
%Test pour la PDMF
%--------------
figure(40)
hold on
[Y, edges] = histcounts(m_tot_pdmf, 500, 'BinLimits',[0,5]);
M = length(Y);
plot(edges(sort([1:M 2:M+1])), Y(sort([1:M 1:M]))/M)
title("PDMF")
set(gca, 'YScale', 'log')
% axis([0 10 1e-2 1e3])

disp(['fraction de rémanents en nombre (WD, NS, BH) ', num2str(mean(frac_N_tot))])
disp(['fraction de rémanents en masse (WD, NS, BH) ', num2str(mean(frac_M_tot))])

figure(41)
[Y_imf, edges] = histcounts(m_tot_imf);
M = length(Y_imf);
plot(edges(sort([1:M 2:M+1])), Y_imf(sort([1:M 1:M]))/M)
title('IMF')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
axis([0 10 1 1e2])

%---------------
%calcul de gamma
%---------------

gamma=tau/uT*2/pi/mean(te)*1e6*365.25;
disp(['gamma (calcule par le te moyen) =    ' num2str(gamma)]);


gam1=4*sqrt(GMsol)/c*uT/sqrt(pc*pc*pc)*length(te)/(n*nbsimul)*86400*365.25*1e6;
gam=gam1*Gammax;
disp(['gamma (integre par MC) = ' num2str(gam)]);


ttobs=tau/(gam/1e6/365.25);
disp(['<tobs> (en jours) = ' num2str(ttobs)]);

taur=gam*pi/2*uT*mean(te)/365.25/1e6;
taur=real(taur);
disp(['tau (avec gamma integré par MC) = ' num2str(taur)]);

%------------------------
% Application du blending
%------------------------

% Param�tre pour le blending 
% f = 0.05;   %fraction des évenements unblendé f = P(1)
% nbar = 4.51; % P(n) = fonction(nbar) = f avec P(n) la proba d'avoir n étoiles dans DeltaS

f = 0.5;
nbar = 1.257;
% f = 0.2;
% nbar = 2.6;
% f = 1;
% nbra = 0;
%retourn teblend (histogramme corrigé) et taurblend (profondeur optique corrigée)
script_blending 

%-----------
%Choix de l'expérience à analyser
%donne teff : tet des observations et eff : efficacité
%------------


temax = 100;
nbre_bin = temax;

%----------------
% Données de l'expérience OGLE
%----------------------

%Exposure
exposure = 2011*table7.N_stars(table7.field == field) /365.25;
% exposure = 2741*sum(table5.N18(extractBetween(table5.field,1,6) == field))*1e-6 /365.25;

%Récupération des évènements du champs

disp(['exposure ogle = ', num2str(exposure)])
disp(['gamma observé ogle = ' num2str(table7.gam(find(extractBetween(table7.field, 1, 6) == field)))])
disp(['tau observé ogle = ' num2str(table7.tau(find(extractBetween(table7.field, 1, 6) == field)))])


%-------------
%Calcul de l'efficacité à partir des données de OGLE IV pour l'article de
%2017 (nature)
%-----------------
figure(1)
eff_2017 = eval(['eff_field.',convertStringsToChars(field)]);

te_inter = 10.^eff_field.log_tE;

te_inter_plot = 10.^[-1; eff_field.log_tE + 0.07];
M = length(eff_2017);
loglog([te_inter_plot(sort([1:M 1:M]));10^(2.5)], [eff_2017(sort([1:M 1:M 1]))])
hold on

teffmaxm=max(te_inter);
teffminm=min(te_inter);

i1_unblend = find((te<=teffmaxm)&(te>=teffminm));
i1_blend = find((teblend<=teffmaxm)&(teblend>=teffminm));

eff_unblend = zeros(1,length(te));	% applique une efficacite nulle aux durees superieures et inferieures
eff_blend = zeros(1,length(teblend));

eff_unblend(i1_unblend) = interp1(te_inter,eff_2017,te(i1_unblend));
eff_blend(i1_blend) = interp1(te_inter,eff_2017,teblend(i1_blend));


%---------------
%efficacité du dossier /eff
%------------------------
VarNames_eff_IV = {'log_tE_min', 'log_tE_max', 'efficiency'};
VarType_eff_IV = {'double', 'double', 'double'};

opts_eff = delimitedTextImportOptions('VariableNames',VarNames_eff_IV,'VariableTypes',VarType_eff_IV,...
                            'Delimiter',delimiter, 'DataLines', 5, ...
                   'WhiteSpace', ' ', 'ConsecutiveDelimitersRule', 'join');
               
eff_field_2019 = readtable(strcat("../OGLEIV/eff/", field, ".eff"),opts_eff);


% eff_2019 = eff_field_2019.efficiency;

%Graph
M = length(eff_field_2019.log_tE_min);
te_graph = [eff_field_2019.log_tE_min(1) ; eff_field_2019.log_tE_min ; eff_field_2019.log_tE_max; eff_field_2019.log_tE_max(M)];
eff_graph = [0 ; eff_field_2019.efficiency(sort([1:M 1:M])) ; 0];
loglog(10.^sort(te_graph), eff_graph*0.4)
legend('2017', '2019 (x0.4)')
title(strcat('efficiency according to the 2 data of the field : ', field));


te_inter_min = 10.^(eff_field_2019.log_tE_min);
te_inter_max = 10.^(eff_field_2019.log_tE_max);

teffmaxm=max(te_inter_max);
teffminm=min(te_inter_min);

% eff_unblend = zeros(1,length(te));	% applique une efficacite nulle aux durees superieures et inferieures
% eff_blend = zeros(1,length(teblend));
% 
% for i = 1:length(te_inter_min)
%     i1_unblend = find(te>=te_inter_min(i) & te<=te_inter_max(i));
%     i1_blend = find(teblend>=te_inter_min(i) & teblend<=te_inter_max(i));
%     
%     eff_unblend(i1_unblend) = ones(size(i1_unblend)) .* eff_field_2019.efficiency(i);
%     eff_blend(i1_blend) = ones(size(i1_blend)) .* eff_field_2019.efficiency(i);
%     
% end


%--------------------------------------------------------------------------------------------------------------------------
% compare le nombre aleatoire precedent a l'efficacite que l'on vient de calculer afin de decider si l'evt est garde ou non
%--------------------------------------------------------------------------------------------------------------------------

%choix de l'éfficacité
eff = eff_2017;

%tirage au sort pour l'efficacité

ra_unblend = rand(1,length(te))*max(eff);
ra_blend = rand(1,length(teblend))*max(eff);

% On choisit l'efficacité ici en prenant les bons indices i

i = find(ra_unblend-eff_unblend<=0); 
teobs = te(i);

ib = find(ra_blend-eff_blend<=0); 
teobsblend = teblend(ib); % On récupère les éléments qui sont soumis au blending avec le calcul d'avant


%---------------
%calcul de gamma
%---------------

disp('Grandeurs avec intervention de l''efficacite experimentale :')

%------------------------------------------------------------------------------------------------
% on ne peut pas calculer le gamma par la formule avec le tobs, car le tau ne prend pas en compte
% l'efficacite. Par contre, on peut deduire tau experimental a partir du gamma calcule par MC
%------------------------------------------------------------------------------------------------

gamobs = gam/length(te)*length(teobs);
disp(['gamma (integre par MC) = ' num2str(gamobs)]);
% 
tauobs=gamobs*pi/2*uT*mean(teobs)/365.25/1e6;
disp(['tau obs (calcule par le te moyen) = ' num2str(tauobs)]);
% 
% tauobsblend=tauobs * gmean * (nbar/(1-exp(-nbar)));
% tauobsblend=real(tauobsblend);
% % disp(['tau observé avec blending (Alibert 2005)  = ' num2str(tauobsblend)]);
% 
gamobsb = gam/length(te)*length(teobsblend);
disp(['gamma avec blending (integre par MC) = ' num2str(gamobsb)]);
% 
tauobsb=gamobsb*pi/2*uT*mean(teobsblend)/365.25/1e6;
disp(['tau obs avec blending (calcule par le te moyen) = ' num2str(tauobsb)]);
% 
% disp(['rapport tau_blend/tau_obs_théorique = ' num2str(tauobsb/tauobs)]);

%------------------------
% affichage des resultats
%------------------------


%telechargement de la courbe du modèle

load ../graph/evenements_1.txt
te_model = evenements_1(:,5);

%Paramètre graph
bin_max = 100;
nbre_bin = bin_max/2;

%Trace la distrib de te  pour le modèle et la courbe stockée localement
[hist, edges] = histcounts(te, nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');
[hist_model, edges] = histcounts(te_model, nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');

%tracé distribution avec blending uniquement la courbe stockée localement
[histb, edges] = histcounts(teblend, nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');

%Courbe expérimentale (avec l'efficacité) :
[hist_obs, edges] = histcounts(teobs, nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');
[hist_obs_b, edges] = histcounts(teobsblend, nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');

centre = zeros(size(edges)-[0,1]);

for j =1:length(centre);
centre(j)=(edges(j)+edges(j+1))/2;
end

%expérience
[hist_exp_normalise, edges] = histcounts(teff, nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');
[hist_exp, edges] = histcounts(teff, nbre_bin, 'BinLimits',[0,bin_max]);


%Graph normalisé
figure(16)
hold on;
plot(centre, hist, 'black');
plot(centre, hist_model, 'red');
legend('local', 'model')
title('comparaison local et modèle')
xlabel('t_{e}')
ylabel('Nombre d''évènements par unité de t_{e}')

%graph en fonction de l'exposition
figure(17)
hold on;
plot(centre, hist_obs.*gamobs*exposure, 'red');
plot(centre, hist_obs_b*gamobsb*exposure, 'black');
plot(centre, hist_exp)
legend('hist modèle', 'hist modèle avec blending (f=0.5)', strcat('OGLE IV,  ', field))
xlabel('t_{e}')
ylabel('Nombre d''évènements par unité de t_{e}')

%Graph normalisé expérience et exp simulée avec l'efficacité pour OGLE IV
figure(18)
hold on;
plot(centre, hist_obs, 'red');
plot(centre, hist_obs_b, 'black');
plot(centre, hist_exp_normalise)
legend('hist modèle', 'hist modèle avec blending (f=0.5)', strcat('OGLE IV,  ', field))
xlabel('t_{e}')
ylabel('Nombre d''évènements par unité de t_{e}')


%%
%Graph normalisé expérience et exp simulée avec l'efficacité pour OGLE III
% figure(18)
% hold on;
% plot(centre, hist_obs, 'red');
% plot(centre, hist_obs_b, 'black');
% M = length(hist_exp_err);
% plot(edges(sort([1:M 1:M])), [0 , 0, hist_exp_err(sort([1:M 2:M-1]))])
% M = length(hist_exp_BW);
% plot(edges(sort([1:M 1:M])), [0 , 0, hist_exp_BW(sort([1:M 2:M-1]))], 'g')
% legend('hist modèle', 'hist modèle avec blending (f=0.5)', 'OGLE III (all stars)', 'OGLE III (fenêtre de Baade)')
% xlabel('t_{e}')
% ylabel('Nombre d''évènements par unité de t_{e}')

%graph noramlisé avec blending
% figure(1);
% hold on;
% plot(centre, hist, 'red')
% plot(centre, histb, 'black')
% title('Blending black et sans blending rouge)');
% hold off;
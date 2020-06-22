%But calculer tau, gamma et <te> pour chaque champ de OGLEIV selon un modèle considéré
%On charge les données OGLE IV et ensuite on calcule pour chaque direction (=champ) via MC ces données

clear

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


fichier ='../graph/comp_modele_OGLEIV.txt';

fid = fopen(fichier, 'w')

fprintf(fid, 'Calcul pour chaque champ des données de mon modèle, comparée avec celle de OGLE IV \n
Les variables stockées sont, dans l'ordre : \n 
Champ; \n longitude; \n lattitude; \n N_stars; \n N_events; 
\n gam_OGLE; \n tau_OGLE; \n mean_te_ogle; \n 
gamma calculé par le modèle (brut); \n tau ----------------- \n mean_te ----------------------\n 
gamma avec efficacité du champ considéré \n tau ---------------- \n mean_te -------------------\n
gamma avec efficacité et blending \n tau ---------------- \n mean_te --------------------')


% for field = table7.field

for field = ["BLG500"; "BLG513"]
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

l = table7.glon(table7.field == field) *pi/180;    % direction d'observation en radian
b = table7.glat(table7.field == field) *pi/180;

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


%-------------------------------------------------------
%calcul de la masse moyenne et de la normalisation de fm
%-------------------------------------------------------

    %--------------------
    % cas du disque mince
    %--------------------

global normfmdm mmeandm

normfmdm=integral(@fmdm,minfdm,msupdm);
mmeandm=integral(@mPmdm,minfdm,msupdm);

    %--------------------
    % cas du disque epais
    %--------------------

global normfmde mmeande

normfmde=integral(@fmde,minfde,msupde);
mmeande=integral(@mPmde,minfde,msupde);


    %-------------
    % cas du bulbe
    %-------------

global normfmbu mmeanbu

normfmbu=integral(@fmbu,minfbu,msupbu);
mmeanbu=integral(@mPmbu,minfbu,msupbu);

    %------------
    % cas du halo
    %------------

    global normfmh mmeanh

normfmh=integral(@fmh,minfh,msuph);
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

% randn('seed',sum(100000*clock)), rand('seed',sum(100000*clock))

x=rand(1,n);
ds=rand(1,n)*(dsup-dinf)+dinf;
m=rand(1,n)*(msuptot-minftot)+minftot;
v=rand(1,n)*vlimit;

gmax(compteur)=max(dgam(x,ds,v,m));

clear x ds m v

end

Gammax=max(gmax);


%---------------------
%---------------------
% debut du Monte-Carlo
%---------------------
%---------------------


%----------------------------------
% ouverture des fichiers de donnees
%----------------------------------

dgamma=[];
dgammaccepte=[];        % initialisation pour les traitement ulterieur des donnees
tecorrespondant=[];

for compteur = 1:nbsimul,

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

%---------------------------------
%tirage de la masse de la lentille
%---------------------------------

ra=rand(1,n);
m(idml) = interp1(ifmmdm,mmdm,ra(idml));
m(idel) = interp1(ifmmde,mmde,ra(idel));
m(ibul) = interp1(ifmmbu,mmbu,ra(ibul));
m(ihl) = interp1(ifmmh,mmh,ra(ihl));

%-------------------------------------------------------
%cas particulier : population de WD dans le de
%-------------------------------------------------------

fraction=0;
raWD=rand(size(idel));
iWD=find(raWD<fraction);
if (length(iWD)>=1)
m(idel(iWD))=0.6*ones(size(iWD));
end;  

%-------------------------------------------------------
%cas particulier : population de WD dans le bulbe
%-------------------------------------------------------

fraction=0;
raWD=rand(size(ibul));
iWD=find(raWD<fraction);
if (length(iWD)>=1)
m(ibul(iWD))=0.6*ones(size(iWD));
end;  

%-------------------------------------------------------
%cas particulier : population de WD dans le dm
%-------------------------------------------------------

fraction=0.;
raWD=rand(size(idml));
iWD=find(raWD<fraction);
if (length(iWD)>=1)
m(idml(iWD))=0.6*ones(size(iWD));
end;  


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
end;


clear x glr glt glp sigrl sigtl sigpl vrotl ds gsr gst gsp sigrs sigts sigps vrots 
clear v m ibul idel idml ihl ibus ides idms ihs ra R z th rhotot evnts te strange

end



%---------------------
%% Analyse des résultats
%---------------------

close all;

%---------------
%calcul de gamma
%---------------

gamma=tau/uT*2/pi/mean(te)*1e6*365.25;
% disp(['gamma (calcule par le te moyen) =    ' num2str(gamma)]);


gam1=4*sqrt(GMsol)/c*uT/sqrt(pc*pc*pc)*length(te)/(n*nbsimul)*86400*365.25*1e6;
gam=gam1*Gammax;
% disp(['gamma (integre par MC) = ' num2str(gam)]);


ttobs=tau/(gam/1e6/365.25);
% disp(['<tobs> (en jours) = ' num2str(ttobs)]);

taur=gam*pi/2*uT*mean(te)/365.25/1e6;
taur=real(taur);
% disp(['tau (avec gamma integré par MC) = ' num2str(taur)]);

%------------------------
% Application du blending
%------------------------

% Param�tre pour le blending 
f = 0.05;   %fraction des évenements unblendé f = P(1)
nbar = 4.51; % P(n) = fonction(nbar) = f avec P(n) la proba d'avoir n étoiles dans DeltaS

f = 0.5;
nbar = 1.257;
% f = 0.2;
% nbar = 2.6;
% f = 1;
% nbra = 0;
%retourn teblend (histogramme corrigé) et taurblend (profondeur optique corrigée)
script_blending 


%Calcul de l'efficacité
VarNames_eff_IV = {'log_tE_min', 'log_tE_max', 'efficiency'};
VarType_eff_IV = {'double', 'double', 'double'};

opts_eff = delimitedTextImportOptions('VariableNames',VarNames_eff_IV,'VariableTypes',VarType_eff_IV,...
                            'Delimiter',delimiter, 'DataLines', 5, ...
                   'WhiteSpace', ' ', 'ConsecutiveDelimitersRule', 'join');
               
eff_field = readtable(strcat("../OGLEIV/eff/", field, ".eff"),opts_eff);


eff = eff_field.efficiency;

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

%tirage au sort pour l'efficacité
ra_unblend = rand(1,length(te))*max(eff_field.efficiency);
ra_blend = rand(1,length(teblend))*max(eff_field.efficiency);

% On choisit l'efficacité ici en prenant les bons indices i

i = find(ra_unblend-eff_unblend<=0); 
teobs = te(i);

ib = find(ra_blend-eff_blend<=0); 
teobsblend = teblend(ib); % On récupère les éléments qui sont soumis au blending avec le calcul d'avant

%-------------------------------
%Exploitation des résultats avec efficacité expérimentale
%-------------------------------

gamobs = gam/length(te)*length(teobs);
% disp(['gamma (integre par MC) = ' num2str(gamobs)]);

tauobs=gamobs*pi/2*uT*mean(teobs)/365.25/1e6;
% disp(['tau obs (calcule par le te moyen) = ' num2str(tauobs)]);

gamobsb = gam/length(te)*length(teobsblend);
% disp(['gamma avec blending (integre par MC) = ' num2str(gamobsb)]);

tauobsb=gamobsb*pi/2*uT*mean(teobsblend)/365.25/1e6;
% disp(['tau obs avec blending (calcule par le te moyen) = ' num2str(tauobsb)]);

%----------------------------------
%enregistrement dans le fichier txt
%------------------------------------


evnts = [field; l*180/pi; b*180/pi; table7.N_stars(table7.field == field); table7.N_events(table7.field == field) ; ...
table7.gam(table7.field == field); table7.tau(table7.field == field)); table7.t_E_mean(table7.field == field)); gam; tau; mean(te); gamobs; tauobs; mean(teobs); ...
gamobsb; tauobsb; mean(teblend)];
fprintf(fid,'%12.8f  %12.8f  %12.8f  %12.8f  %12.8f %12.8f  %12.8f  %12.8f %12.8f  %12.8f  %12.8f %12.8f  %12.8f  %12.8f %12.8f  %12.8f  %12.8f\n',evnts);
end
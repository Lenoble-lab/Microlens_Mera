%fonction main_microlens qui retourne dans un fichier les paramètres des évè sélectionné.
%Peut être utilisé directement dans un prgm pour tracer un histogramme
%Utilise le modèle directement mais attention à tous les paramètres

% function [tau]  = PDMF_Maraston_1(vlimit, dsup, dinf, n, nbsimul, Ro, elev, Rcoro, l, b,...
    % uT, AT, minf, msup, Vinf, Vsup)

global vsr vsp vst vlp vlt vlr
%----------------------------------
% Constantes physiques (unites SI)
%----------------------------------
G=6.672e-11;	pc=3.08567802e16;
kpc=pc*1e3;  	Msol=1.989e30;
c=299792458;	GMsol=1.32712497e20;

%----------------
%liste pour observer ensuite les distributions
%----------------------
m_tot_imf = [];
m_tot_pdmf = [];
frac_M_tot = [];
frac_N_tot = [];
v_tot = [];
L_tot = [];
eve_BH = [];
star_pop = [];
frac_eve_tot = [];
%-----------------------------------
% Param�tres de la fonction de masse
%-----------------------------------
global minfbu msupbu minfdm msupdm minfde msupde minfh msuph 

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

global  normfl

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
fds = nsource (dd);
ifds = cumsum(fds);	% primitive 
ifds = ifds-ifds(1);
ifds = real(ifds./ifds(end));	% on fait en sorte que la primitive varie de 0 a 1


%-------------------------------------------
% calcul du maximum de Gamma par monte carlo
%-------------------------------------------

couple_max = fminsearch(@dgammax,[0.5, 15e3]);
Gammax = -dgammax(couple_max);

disp(['couple max : (x,L) = ' num2str(couple_max)])
disp(['Gammax = ' num2str(Gammax*1e-12) '*1e12']);

%%
%---------------------
%---------------------
% debut du Monte-Carlo
%---------------------
%---------------------


%----------------------------------
% ouverture des fichiers de donnees
%----------------------------------

fid1 = fopen('evenements.txt','w');
fids = fopen('strange.txt','w');

dgamma=[];
dgammaccepte=[];        % initialisation pour les traitement ulterieur des donnees
tecorrespondant=[];

for compteur = 1:nbsimul

  if ( rem(compteur,100) == 0 )  
    
disp(['compteur = ' num2str(compteur)]);

  end

%----------------------------------
%tirage de la distance de la source
%----------------------------------

ra = rand(1,n);
[ifds, index] = unique(ifds); 

ds = interp1(ifds,dd(index),ra);
L_tot = [L_tot, ds];

%------------
% Tirage de x
%------------

% dl=rand(1,n).*ds;
% x = dl./ds;
% rapport = dinf./ds;
% x=(rand(1,n)-rapport)/(1-rapport) + rapport;
x = rand(1,n);
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
[ifmmdm, index] = unique(ifmmdm, 'last'); 
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
%star_pop : code [1,2,3,4,5] = [BD, MS, WD, NS, BH]
star_pop = zeros(size(x));

% m_tot_imf = [m_tot_imf, m];
% [m, frac_N, frac_M] = PDMF_gould_1(m);
[m, frac_N, frac_M, frac_eve, star_pop] = PDMF_Maraston_1(m);
frac_N_tot = [frac_N_tot ; frac_N];
frac_M_tot = [frac_M_tot ; frac_M];
frac_eve_tot = [frac_eve_tot; frac_eve];

% m = ones(size(x));
m_tot_pdmf = [m_tot_pdmf, m];

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
% % sigps(ihs)=sigph(R(ihs),z(ihs),th(ihs));



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

v_tot = [v_tot, v];
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
i = find(g>=ra);

j = find(g>1);

dgamma = [dgamma,g];

dgammaccepte=[dgammaccepte,g(i)];

%---------------
%sélection des eve de BH pour test
%--------------------
i_BH_eve = find(star_pop == 5);
te_BH = (2/c*sqrt(GMsol*pc)/86400).*sqrt(ds(i_BH_eve).*m(i_BH_eve).*x(i_BH_eve).*(1-x(i_BH_eve)))./v(i_BH_eve);

eve_BH = [eve_BH, [x(i_BH_eve);ds(i_BH_eve); v(i_BH_eve); m(i_BH_eve); te_BH; star_pop(i_BH_eve); dgam(x(i_BH_eve),ds(i_BH_eve),v(i_BH_eve),m(i_BH_eve))./Gammax]];



%------------------------
%selection des evenements
%------------------------

if (length(i)>=1)
  te = (2/c*sqrt(GMsol*pc)/86400).*sqrt(ds(i).*m(i).*x(i).*(1-x(i)))./v(i);
  
  tecorrespondant=[tecorrespondant,te];
  
  evnts = [x(i);ds(i);v(i);m(i);te;star_pop(i)];
  fprintf(fid1,'%12.8f  %12.8f  %12.8f  %12.8f  %12.8f  %12.8f\n',evnts);
end;


if (length(j)>=1)
  strange = [x(j);ds(j);v(j);m(j)];
  fprintf(fids,'%12.8f  %12.8f  %12.8f  %12.8f\n',strange);
end;


clear x glr glt glp sigrl sigtl sigpl vrotl ds gsr gst gsp sigrs sigts sigps vrots 
clear v m ibul idel idml ihl ibus ides idms ihs ra R z th rhotot evnts te strange star_pop
end
%---------------------
%% Analyse des résultats
%---------------------

fid=fopen('simul_para.txt','w');

fprintf(fid,'%15.11f  \n',tau);

fprintf(fid,'%15.11f  \n',n);
fprintf(fid,'%15.11f  \n',nbsimul);
fprintf(fid,'%15.11f  \n',tau);
fprintf(fid,'%15.11f  \n',Gammax);
fprintf(fid,'%15.11f  \n',uT);
fprintf(fid,'%15.11f  \n',AT);
fprintf(fid, '%15.11f  \n',mean(m_tot_pdmf));

fclose(fid);

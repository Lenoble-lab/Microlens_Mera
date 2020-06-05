% Programme de micro-lentilles gravitationnelles (bulbe et disque)

clear

%------------------
%fichiers resultats
%------------------

fichevents='evenements.sim1.txt';
fichres='resultat.sim1.txt';
fichpara='para.sim1.txt';

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
% n = 10000;
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

% definition de la fenetre de Baade dans la majeure partie des articles :  l = 1 et b = -4
% A priori c'est cette definition qui est juste.
l = 1 *pi/180;    % direction d'observation en radian
b = -4 *pi/180;


% definition de la fenetre de Baade dans les theses de Mera et Alibert : l = 4 et b = -1
l = 4 *pi/180;    % direction d'observation en radian
b = -1 *pi/180;

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
%% Preparation de la table pour le tirage aleatoire de la masse
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
% test : toutes les lentilles ont la meme masse 
%----------------------------------------------

%m=0.6*ones(size(x));

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
fid=fopen(fichres,'w');

fprintf(fid,'%15.11f  \n',tau);
fprintf(fid,'%12.8f  \n',Gammax);

fclose(fid);

fid=fopen(fichpara,'w');

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


%----------------------------------------
%recuperation des evenements selectionnes
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
% N=gam*exposure;
% disp(['nb d''evt  = ' num2str(N)]);

taur=gam*pi/2*uT*mean(te)/365.25/1e6;
taur=real(taur);
disp(['tau (avec gamma integré par MC) = ' num2str(taur)]);


f_n_bar = @(n) n.*exp(-n)./(1-exp(-n));

n = 1e2;
n_li = (0:1e-3:1).*20;

%f_li : fraction d'étoiles non blendée
f_li = f_n_bar(n_li);

for k = 1:length(f_li)
    
nbar = n_li (k);
f = f_li(k);

%retourne teblend (histogramme corrigé) et taurblend (profondeur optique corrigée)
script_blending 

%-----------
%Choix de l'expérience à analyser
%donne teff : te des observations et eff : efficacité
%------------

exp_macho_2005

%------------------------------------------------------------------------------------------------
% on ne peut pas calculer le gamma par la formule avec le tobs, car le tau ne prend pas en compte
% l'efficacite. Par contre, on peut deduire tau experimental a partir du gamma calcule par MC
%------------------------------------------------------------------------------------------------

gamobs = gam/length(te)*length(teobs)*max(eff);

tauobs=gamobs*pi/2*uT*mean(teobs)/365.25/1e6;

gamobsb = gam/length(te)*length(find(teobsblend~=0))*max(eff);

tauobsb=gamobsb*pi/2*uT*mean(teobsblend)/365.25/1e6;

rapport(k) = tauobsb/tauobs;
end

figure(1)
hold on;
plot(f_li, rapport)
ylabel('\tau_{obs} / \tau')
xlabel('fraction de sources sans biais de confusion')
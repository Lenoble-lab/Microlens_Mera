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
elev = 0 ;      % elevation au dessus du plan du disque
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
% l = 1 *pi/180;    % direction d'observation en radian
% b = -4 *pi/180;


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

%----------------------------------
% comparaison avec KP
% toutes les masses sont egales a 0.6
%----------------------------------

%mmeanh=0.6;
%mmeanbu=0.6;
%mmeande=0.6;
%mmeandm=0.6;


disp(['masse moyenne du disque mince = ' num2str(mmeandm)]);
disp(['masse moyenne du disque epais = ' num2str(mmeande)]);
disp(['masse moyenne du bulbe = ' num2str(mmeanbu)]);
disp(['masse moyenne du halo = ' num2str(mmeanh)]);


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

% ra=rand(1,n);
% [R,z,th]= toGC(ds);
% rhotot = rhodm(R,z,th) + rhode(R,z,th) + rhobulbe(R,z,th) + rhohalo(R,z,th);
% idms=find(ra<=(rhodm(R,z,th)./rhotot));
% ides=find(rhodm(R,z,th)./rhotot<ra & ra <= (rhodm(R,z,th)+rhode(R,z,th))./rhotot);
% ibus=find((rhodm(R,z,th)+rhode(R,z,th))./rhotot<ra & ra <= (rhodm(R,z,th)+rhode(R,z,th)+rhobulbe(R,z,th))./rhotot);
% ihs=find(ra >= (rhodm(R,z,th)+rhode(R,z,th)+rhobulbe(R,z,th))./rhotot);

%--------------------------------------------
%a quelle population appartient la lentille ?
%--------------------------------------------

% ra=rand(1,n);
% [R,z,th]= toGC(x.*ds);
% rhotot = rhodm(R,z,th) + rhode(R,z,th) + rhobulbe(R,z,th) + rhohalo(R,z,th);
% idml=find(ra<=(rhodm(R,z,th)./rhotot));
% idel=find(rhodm(R,z,th)./rhotot<ra & ra <= (rhodm(R,z,th)+rhode(R,z,th))./rhotot);
% ibul=find((rhodm(R,z,th)+rhode(R,z,th))./rhotot<ra & ra <= (rhodm(R,z,th)+rhode(R,z,th)+rhobulbe(R,z,th))./rhotot);
% ihl=find(ra >= (rhodm(R,z,th)+rhode(R,z,th)+rhobulbe(R,z,th))./rhotot);


%-----------------------------------------------------------------
% test : methode de Peale. Lentilles et sources ont les proprietes
% du disque si R>Rcoro, et du bulbe sinon
%-----------------------------------------------------------------

[R,z,th]= toGC(x.*ds);
ibul=find(R<=Rcoro);
idml=find(R>Rcoro);
idel=find(R<0.);
ihl=find(R<0.);

[R,z,th]= toGC(ds);
ibus=find(R<=Rcoro);
idms=find(R>Rcoro);
ides=find(R<0.);
ihs=find(R<0.);


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
% sigrl(idel)=sigrde(R(idel),z(idel),th(idel));
% sigtl(idel)=sigtde(R(idel),z(idel),th(idel));

%Bulbe
sigrl(ibul)=sigrb(R(ibul),z(ibul),th(ibul));
sigtl(ibul)=sigtb(R(ibul),z(ibul),th(ibul));

%Halo
% sigrl(ihl)=sigrh(R(ihl),z(ihl),th(ihl));
% sigtl(ihl)=sigth(R(ihl),z(ihl),th(ihl));


%Coordonnée cylindrique u_z
% sigzl(idml)=sigzdm(R(idml),z(idml),th(idml));
% sigzl(ibul)=sigzb(R(ibul),z(ibul),th(ibul));
% % sigzl(idel)=sigzde(R(idel),z(idel),th(idel));
% sigzl(ihl)=sigzh(R(ihl),z(ihl),th(ihl));

%Coordonnée sphérique u_phi
sigpl(idml)=sigpdm(R(idml),z(idml),th(idml));
sigpl(idel)=sigpde(R(idel),z(idel),th(idel));
sigpl(ibul)=sigpb(R(ibul),z(ibul),th(ibul));
sigpl(ihl)=sigph(R(ihl),z(ihl),th(ihl));


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
% sigrs(ides)=sigrde(R(ides),z(ides),th(ides));
% sigts(ides)=sigtde(R(ides),z(ides),th(ides));

%Bulbe
sigrs(ibus)=sigrb(R(ibus),z(ibus),th(ibus));
sigts(ibus)=sigtb(R(ibus),z(ibus),th(ibus));

%Halo
% sigrs(ihs)=sigrh(R(ihs),z(ihs),th(ihs));
% sigts(ihs)=sigth(R(ihs),z(ihs),th(ihs));

%Coordonnée cylindrique u_z
% sigzs(idms)=sigzdm(R(idms),z(idms),th(idms));
% sigzs(ides)=sigzde(R(ides),z(ides),th(ides));
% sigzs(ibus)=sigzb(R(ibus),z(ibus),th(ibus));
% sigzs(ihs)=sigzh(R(ihs),z(ihs),th(ihs));


%Coordonnée sphérique u_phi
sigps(idms)=sigpdm(R(idms),z(idms),th(idms));
sigps(ides)=sigpde(R(ides),z(ides),th(ides));
sigps(ibus)=sigpb(R(ibus),z(ibus),th(ibus));
sigps(ihs)=sigph(R(ihs),z(ihs),th(ihs));



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
v = vperp(x,glr,glt,glp,sigrl,sigtl,sigpl,vrotl,ds,gsr,gst,gsp,sigrs,sigts,sigps,vrots);


%coordonées cylindriques
% v_cyl = vperp_cyl(x,glr,glp,glt,sigrl,sigtl,sigzl,vrotl,ds,gsr,gsp,gst,sigrs,sigts,sigzs,vrots);

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



%if ( compteur == nbsimul )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TEST des distributions des VA principales

%figure(1);
%hist(x,100);
%mx=mean(x);
%sx=std(x);
%title('Distribution de la VA x pendant le Monte-carlo');
%hold off;
%OK : loi uniforme
%figure(2);
%hist(ds,100);
%mds=mean(ds);
%sds=std(ds);
%title('Distribution de la VA ds pendant le Monte-carlo');
%hold off;
% OK donne une distribution correcte max sur le GC
%figure(3);
%hist(v,100);
%mv=mean(v);
%sv=std(v);
%title('Distribution de la VA v pendant le Monte-carlo');
%hold off;
%%distrib correcte???
%figure(4);
%hist(m,10000);
%mm=mean(m);
%sm=std(m);
%title('Distribution de la VA m pendant le Monte-carlo');
%hold off;
% OK: loi de la FM (log-normal et expo)
%end;


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
disp(['gamma (calcule par le te moyen) = ' num2str(gamma)]);


% Juste un ou deux tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ttobs=tau/(gamma/1e6/365.25);
disp(['<tobs1> (en jours) = ' num2str(ttobs)]);
t=pi/2*te;
l=length(t);
l;
disp(['nombre d''evenements 1 = ' num2str(l)]);
ttobs2=sum(t)/l;
disp(['<tobs2> (en jours) = ' num2str(ttobs2)]);
%c'est le meme resultat par definition de mean()
expos=(sum(t)/tau)/(1e6*365.25);
disp(['Exposure (en 10^6 an-etoile) = ' num2str(expos)]);
N=(gamma)*expos;
disp(['nombre d''evenements 2 = ' num2str(N)]);



gam=4*sqrt(GMsol)/c*uT/sqrt(pc*pc*pc)*length(te)/(n*nbsimul)*86400*365.25*1e6*Gammax;
disp(['gamma (integre par MC) = ' num2str(gam)]);



ttobs=tau/(gam/1e6/365.25);
disp(['<tobs> (en jours) = ' num2str(ttobs)]);
N=gam*expos;
disp(['nb d''evt  = ' num2str(N)]);


%------------------------------------
%------------------------------------
% application du facteur d'efficacite
%------------------------------------
%------------------------------------

%------------------------
% Efficacit�s MACHO bulbe
%------------------------



% Durees observees en direction du bulbe, incluant les 4 evts longs
tebulbe = [4.4000,5.5000,5.5500,6.2000,6.3500,6.7000,7.2000,7.4000,7.5000,7.9000,8.2000,8.2000,8.2500,10.9500,10.9500,10.9500,10.9500,11.7500,11.8500,12.3000,14.7500,14.8000,14.9500,15.3500,16.0500,16.2000,16.4500,16.9500,18.1000,20.5000,21.3000,22.2500,24.2500,24.7000,26.9500,27.2000,27.7500,33.3500,75.8500,75.9000,97.8500,112.2500];
% Pour ce tebulbe les donnees sont deja des rayons d'einstein
sampleffA=[0.0982,0.1240,0.1266,0.1395,0.1434,0.1498,0.16015,0.1643,0.1653,0.1689,0.1741,0.1741,0.1746,0.2015,0.2015,0.2015,0.2015,0.2113,0.2144,0.2170,0.2390,0.2392,0.2402,0.2449,0.2526,0.2537,0.2547,0.2583,0.2645,0.2790,0.2841,0.2893,0.3017,0.3048,0.3141,0.3159,0.3177,0.3286,0.4001,0.4006,0.4102,0.4164];
sampleffB=[0.0891,0.1214,0.1240,0.14465,0.1498,0.1576,0.1710,0.1782,0.1818,0.1834,0.1932,0.1932,0.1937,0.2376,0.2376,0.2376,0.2376,0.2475,0.2531,0.2557,0.2819,0.2821,0.2854,0.2908,0.2991,0.3007,0.3022,0.3058,0.3125,0.3306,0.3384,0.3430,0.3565,0.3601,0.3694,0.3720,0.3745,0.3880,0.4450,0.4456,0.4489,0.4494];
photomeffA=[0.0917,0.1111,0.1147,0.1276,0.1307,0.1343,0.1472,0.1498,0.1500,0.1576,0.1632,0.1632,0.1640,0.2066,0.2066,0.2066,0.2066,0.2170,0.2221,0.2273,0.2542,0.2547,0.2552,0.2609,0.2645,0.2666,0.2686,0.2717,0.2769,0.2919,0.2965,0.3007,0.3100,0.3151,0.3220,0.3255,0.3280,0.3523,0.4190,0.4196,0.4262,0.4288];
photomeffB=[0.0826,0.1111,0.1147,0.1317,0.1421,0.1493,0.1679,0.1756,0.1780,0.1860,0.1937,0.1937,0.1950,0.2635,0.2635,0.2635,0.2635,0.2790,0.2841,0.2945,0.3332,0.3345,0.3358,0.3410,0.3461,0.3490,0.3513,0.3580,0.3668,0.3875,0.3968,0.40295,0.4185,0.4288,0.4350,0.4443,0.4520,0.4644,0.5364,0.5370,0.5321,0.5350];


tmachob = [ 0, 0.56, 0.7049, 0.8872, 1.116, 1.405, 1.769, 2.227, 2.803, 3.540, 4.456, 5.609, 7.060, 8.887, 11.18, 14.08, 17.72, 22.3, 28.08, 35.34, 44.49, 56, 70.49, 89.02, 112.05, 141.04];
tmachob = tmachob/2.; %ATTENTION, lequel est le bon?????
stdeffmachob = [ 0, 0, 0, 0.00091407, 0.016453, 0.0338, 0.05758, 0.08135, 0.11791, 0.14899, 0.1718, 0.2065, 0.2404, 0.2815, 0.3071, 0.3436, 0.3839, 0.4076, 0.4095, 0.4369, 0.4716, 0.4561, 0.4570, 0.4360, 0.4140, 0.3884];
clpeffmachob = [ 0, 0.0054844, 0.001828, 0.01005, 0.08866, 0.1809, 0.2568, 0.3171, 0.4085, 0.4515, 0.5018, 0.5237, 0.5840, 0.6087, 0.6590, 0.6636, 0.7129, 0.7248, 0.8089, 0.7385, 0.7888, 0.7641, 0.6974, 0.6663, 0.6069, 0.5475];
% ces quatres dernier tableaux sont tires des figures 1 et 12 de l'article The Macho Project :  
sampleffmachobA = [0,0,0,0,0,0.002,0.0052,0.0103,0.0188,0.0284,0.0403,0.0589,0.0646,0.0896,0.1240,0.1563,0.1808,0.2041,0.2299,0.2609,0.2893,0.3182,0.3317,0.3616,0.3823,0.3961];
sampleffmachobB = [0,0,0,0,0,0.001,0.0013,0.0026,0.0065,0.0129,0.0222,0.0413,0.0788,0.0992,0.1214,0.1679,0.2015,0.2392,0.2712,0.3100,0.3435,0.3771,0.3926,0.4197,0.4365,0.4435];
photomeffmachobA = [0,0,0,0,0,0.0008,0.002,0.0057,0.0103,0.0181,0.0315,0.0459,0.0671,0.0916,0.1151,0.1471,0.1755,0.2104,0.2452,0.2757,0.3004,0.3317,0.3614,0.3872,0.4027,0.4181];
photomeffmachobB = [0,0,0,0,0,0,0.001,0.0021,0.0044,0.0090,0.0155,0.0399,0.0542,0.0839,0.1151,0.1673,0.2117,0.2684,0.3175,0.3614,0.3985,0.4517,0.4698,0.5033,0.5255,0.5384];



%-------------------------------
% Efficacites MACHO 24 Mars 2000
% Criteres A ou B
%-------------------------------


%attention, les machos prennent le diametre d'einstein --> diviser les tmacho par 2 !


tmachoNEW = [ 1, 1.46, 2, 2.7, 3.8, 5.2, 7.1, 9.7, 13.3, 18.2, 24.9, 34.2, 46.8, 64.2, 88.0, 121.0, 165.0, 227.0, 311.0, 426.0, 584.0, 800.0, 1096.0, 1502.0 ];
tmachoNEW = tmachoNEW/2.;
effmachoA = [8.7e-4, 2.1e-3, 6.5e-3, 1.3e-2, 3.2e-2, 5.2e-2, 8.7e-2, 1.1e-1, 1.5e-1, 1.9e-1, 2.5e-1, 2.8e-1, 3.2e-1, 3.6e-1, 3.9e-1, 4.1e-1, 4.2e-1, 4.3e-1, 4.1e-1, 3.2e-1, 1.0e-1, 4.7e-2, 3.3e-2, 3.0e-2 ];
effmachoB = [1.9e-4, 5.4e-4, 2.2e-3, 5.7e-3, 1.5e-2, 3.7e-2, 7.8e-2, 1.1e-1, 1.8e-1, 2.5e-1, 3.2e-1, 3.8e-1, 4.4e-1, 4.7e-1, 5.1e-1, 5.4e-1, 5.3e-1, 5.3e-1, 4.9e-1, 4.1e-1, 1.7e-1, 9.1e-2, 6.8e-2, 6.1e-2 ]; 

% ces quatres dernier tableaux sont tires des figures 1 et 12 de l'article The Macho Project : microlensing detection efficiency de Alcock et al, date du 21 march 2001.
sampleffmachoA = [0,0.0025,0.0077,0.0176,0.032,0.0532,0.0801,0.1100,0.1472,0.1834,0.2185,0.2561,0.2970,0.3270,0.3565,0.3875,0.4035,0.4174,0.4236,0.3802,0.1756,0.0093,0,0];
sampleffmachoB = [0,0.001,0.002,0.0057,0.0155,0.0362,0.0646,0.1007,0.1550,0.2072,0.2562,0.3027,0.3487,0.3875,0.4159,0.4391,0.4469,0.4494,0.4391,0.3934,0.1808,0.0093,0,0];
photomeffmachoA = [0,0,0.004,0.0095,0.0222,0.0426,0.0681,0.0991,0.1353,0.1760,0.2277,0.2710,0.3072,0.3485,0.3820,0.4078,0.4220,0.4300,0.4187,0.3655,0.2117,0.0718,0.0390,0.0315];
photomeffmachoB= [0,0.0012,0.0023,0.0045,0.0103,0.0258,0.0542,0.0955,0.1523,0.2168,0.2937,0.3562,0.4130,0.4620,0.4945,0.5291,0.5356,0.5322,0.4976,0.4228,0.2839,0.1208,0.0779,0.0645];


%TRACE DES EFFICACITES EN FONCTION DU TEMPS (ECHELLE LOG) 

% lt=log14(tmachob);
% figure(5);
% title('graphe des differentes efficacites pour les donnees machob');
% hold on;
% plot(lt,stdeffmachob,'g--');
% plot(lt,clpeffmachob,'c--');
% plot(lt,sampleffmachobA,'b-');
% plot(lt,sampleffmachobB,'r-');
% plot(lt,photomeffmachobA,'b-.');
% plot(lt,photomeffmachobB,'r-.');
% legend('stdeff','clpeff','sampleffA','sampleffB','photomeffA','photomeffB');
% hold off;
% 
% lt=log10(tebulbe);
% figure(51);
% title('graphe des differentes efficacites pour les donnees tebulbe');
% hold on;
% plot(lt,sampleffA,'b-');
% plot(lt,sampleffB,'r-');
% plot(lt,photomeffA,'b-.');
% plot(lt,photomeffB,'r-.');
% legend('sampleffA','sampleffB','photomeffA','photomeffB');
% hold off;
% 
% lt=log10(tmachoNEW);
% figure(6);
% title('graphe des differentes efficacites pour les donnees machoNEW');
% hold on;
% plot(lt,effmachoA,'g--');
% plot(lt,effmachoB,'c--');
% plot(lt,sampleffmachoA,'b-');
% plot(lt,sampleffmachoB,'r-');
% plot(lt,photomeffmachoA,'b-.');
% plot(lt,photomeffmachoB,'r-.');
% legend('effA','effB','sampleffA','sampleffB','photomeffA','photomeffB');
% hold off;


%----------------------
% Choix de l'efficacite
%----------------------

eff = sampleffA; 

teff = tebulbe;


%----------------------
% essai, efficacite = 1
%----------------------

%eff=ones(size(teff));


%-----------------------------------------------------------------------------------------------
% Interpolation lineaire de l'efficacite pour determiner la probabilite qu'un evt a d'etre garde
%-----------------------------------------------------------------------------------------------

%if ((eff==stdeffmachob)|(eff==clpeffmachob)|(eff==effmachoA)|(eff==effmachoB)|(eff==ones(size(teff))))
%    tinterp=teff;
%    effinterp=eff;
%else
%    tinterp=tmachoNEW; 
%    if ((eff==sampleffA)|(eff==sampleffmachobA)|(eff==sampleffmachoA))
%        effinterp=sampleffmachoA;
%    elseif ((eff==sampleffB)|(eff==sampleffmachobB)|(eff==sampleffmachoB))
%        effinterp=sampleffmachoB;
%    elseif ((eff==photomeffA)|(eff==photomeffmachobA)|(eff==photomeffmachoA))
%        effinterp=photomeffmachoA;
%    elseif ((eff==photomeffB)|(eff==photomeffmachobB)|(eff==photomeffmachoB))
%        effinterp=photomeffmachoB;
%    end;
%end;

tinterp=tmachoNEW;
effinterp=sampleffmachoA;

teffmax=max(tinterp);
teffmin=min(tinterp);

i1 = find((te<=teffmax)&(te>=teffmin));
effsim = zeros(1,length(te));	% applique une efficacite nulle aux durees superieures et inferieures
effsim(i1) = interp1(tinterp,effinterp,te(i1));


%-----------------------------------------------------------------------------
% Tirage d'un nombre aleatoire qui servira a decider si l'evt est garde ou non
%-----------------------------------------------------------------------------

ra = rand(1,length(te))*max(effinterp);


% Trace dgamma accepte en fonction de tecorrespondant

size(tecorrespondant);   % verification de leurs tailles pour le tracage
% size(dgammaccepte);
% figure(7);
% plot(tecorrespondant,dgammaccepte,'b.');
% title('dgammaccepte en fonction de tecorrespondant (resultat du Monte-Carlo)');
% hold off;


ecar=0.5;
for i = ecar:ecar:max(tecorrespondant);
    j=find( (i-ecar)<tecorrespondant & tecorrespondant<=i );
    k=floor(2*i);
    d=dgammaccepte(j);
    d=[d,0]; % Pour qu'il y ait au moins un elements (cas de la matrice vide)
    dg(k)=max(d);
end;
i = ecar:ecar:max(tecorrespondant);
% figure(8);
% plot(i,dg,'r-');
% axis([0,100,0,1]);
% title('dgammaccepte en fonction de tecorrespondant ( envelope ou fonction reelle )');
% hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(9);
% ep=interp1(tinterp,effinterp,i);
% des=ep.*dg;
% plot(i,des,'g-');
% axis([0,100,0,0.05]);
% title('taux d''evt en fonction de te (a un facteur pres  ( epsilon*dgamma) )');
% hold off;



% % Trace de la distrib de te (simulation)
% figure(10);
% hist(te,1000);
% title('Distribution de te obtenu a partir de la simulation');
% hold off;
% 
% %Trace de la fonction d'efficacite
% figure(11);
% lt = log(te);
% plot(lt,effsim,'.');
% title('efficacite de la simulation en fonction des durees d''evenements (echelle log)');
% hold off;
% 
% 

%--------------------------------------------------------------------------------------------------------------------------
% compare le nombre aleatoire precedent a l'efficacite que l'on vient de calculer afin de decider si l'evt est garde ou non
%--------------------------------------------------------------------------------------------------------------------------

i1 = find(ra-effsim<=0);
teobs=te(i1);
% 
% lt=log10(teobs);
% figure(81); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot(lt,ra(i1),'.');
% title('resultat du Monte-Carlo sur l''efficacite');
% hold off;

%teobs=te;  c'est faux puisqu'il faut faire
%l'integration qui se fait soit par le monte carlo soit en essayant de
%faire le veritable calcul (un des autres essai que j'ai fait).


%---------------
%calcul de gamma
%---------------

disp('Grandeurs avec intervention de l''efficacite experimentale :')

%------------------------------------------------------------------------------------------------
% on ne peut pas calculer le gamma par la formule avec le tobs, car le tau ne prend pas en compte
% l'efficacite. Par contre, on peut deduire tau experimental a partir du gamma calcule par MC
%------------------------------------------------------------------------------------------------

gamobs = gam/length(te)*length(teobs)*max(eff);
disp(['gamma (integre par MC) = ' num2str(gamobs)]);

tauobs=gamobs*pi/2*uT*mean(teobs)/365.25/1e6;
disp(['tau obs (calcule par le te moyen) = ' num2str(tauobs)]);


%------------------------
%------------------------
% affichage des resultats
%------------------------
%------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Test de l'influence du MC sur les distributions

%figure(12);
%hist(x,100);
%mx=mean(x);
%sx=std(x);
%title('Distribution de la VA x apres le Monte-carlo (influence) ');
%hold off;
%%EST-CE QU'ON PEUT DIRE QUE C'EST LA REPRESENTATION DE LA LOI dGAMMA???
%figure(13);
%hist(ds,100);
%mds=mean(ds);
%sds=std(ds);
%title('Distribution de la VA ds apres le Monte-carlo (influence) ');
%hold off;
%figure(14);
%hist(v,100);
%mv=mean(v);
%sv=std(v);
%title('Distribution de la VA v apres le Monte-carlo (influence) ');
%hold off;
%figure(15);
%hist(m,100000);
%mm=mean(m);
%sm=std(m);
%title('Distribution de la VA m apres le Monte-carlo (influence) ');
%hold off;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%Trace la distribde te en zoomant
[hist, edges] = histcounts(te, 3'0, 'BinLimits',[0,30], 'Normalization', 'probability');
centre = zeros(size(edges)-[0,1]);

for j =1:length(centre);
centre(j)=(edges(j)+edges(j+1))/2;
end

figure(18);
% histogram(te, 30, 'BinLimits',[1,30]);
plot(centre, hist);
title('Distribution de te (en zoomant)');
hold off;

%trace la distrib de teobs
figure(17);
hist(teobs,1000);
title('Distribution de teobs (evenements garde apres le dernier Monte-Carlo)');
hold off;




temaxgraphe=100;  % pour le trace des graphe
temax=5000;       % pour le calcul integral

% Definition des differents nombres de bin pour les graphiques
nbbin1=temaxgraphe;
nbbin2=temaxgraphe/2;     
nbbin5=temaxgraphe/5;

bords1 = zeros(1, nbbin1) ; centre1 = zeros(1, nbbin1-1);
bords2 = zeros(1, nbbin2) ; centre2 = zeros(1, nbbin2-1);
bords5 = zeros(1, nbbin5) ; centre5 = zeros(1, nbbin5-1);

exposure=6.5244;
teobs = te;
%Definition des tableaux bords et centre pour une largeur de 1 jour
for j =1:nbbin1;
bords1(j)=temaxgraphe/nbbin1*(j-1);
end;
for j =1:nbbin1-1;
centre1(j)=(bords1(j)+bords1(j+1))/2;
end
centre1(nbbin1)=bords1(nbbin1)+(bords1(nbbin1)-bords1(nbbin1-1))/2;

i=find(teobs<temaxgraphe); % ce tableau i est le meme pour les trois graphes

h1=histc(teobs(i),bords1)*gamobs*exposure/length(teobs(i));

%Definition des tableaux bords et centre pour une largeur de 2 jour
for j =1:nbbin2;
bords2(j)=temaxgraphe/nbbin2*(j-1);
end;
for j =1:nbbin2-1;
centre2(j)=(bords2(j)+bords2(j+1))/2;
end
centre2(nbbin2)=bords2(nbbin2)+(bords2(nbbin2)-bords2(nbbin2-1))/2;

h2=histc(teobs(i),bords2)*gamobs*exposure/length(teobs(i));

%Definition des tableaux bords et centre pour une largeur de 5 jour
for j =1:nbbin5;
bords5(j)=temaxgraphe/nbbin5*(j-1);
end;
for j =1:nbbin5-1;
centre5(j)=(bords5(j)+bords5(j+1))/2;
end
centre5(nbbin5)=bords5(nbbin5)+(bords5(nbbin5)-bords5(nbbin5-1))/2;

h5=histc(teobs(i),bords5)*gamobs*exposure/length(teobs(i));


% Tracage de ces trois graphes ensemble
figure(18);
hold off;
hold on;
plot(centre1,h1,'b-');
plot(centre2,h2,'r-');
plot(centre5,h5,'g-');
axis([0,30,0,20]);
title('Diagramme des evenements obtenu par simulation pour des largeurs de representation differentes');
legend('largeur de 1 jour','largeur de 2 jours','largeur de 5 jours');
xlabel('durees d''evenements');
ylabel('nombre d''evenements par unite de temps');
hold off;


% TRACAGE des trois graphiques precedents separement avec la distribution experimentale correspondante

figure(19);
j=find(teff<temaxgraphe);  % c'est le meme tableau pour les trois graphes
hist(teff(j),bords1);% essayer aussi hi1=hist(teff(i),bords1)*gamobs*exposure/length(teff(i));
teff(j);
%bar(centre1,hi1,'r');
hold on;
plot(centre1,h1,'r');
title('comparaison experience/simulation pour une largeur de 1 jours');
legend('experience','simulation');
hold off;

figure(20);
hist(teff(j),bords2);   % essayer aussi hi2=hist(teff(i),bords2)*gamobs*exposure/length(teff(i));
%bar(centre2,hi2,'r');
hold on;
plot(centre2,h2,'r');
title('comparaison experience/simulation pour une largeur de 2 jours');
legend('experience','simulation');
hold off;

figure(21);
hist(teff(j),bords5);   % essayer aussi hi5=hist(teff(i),bords5)*gamobs*exposure/length(teff(i));
%bar(centre5,hi5,'r');
hold on;
plot(centre5,h5,'r');
title('comparaison experience/simulation pour une largeur de 5 jours');
legend('experience','simulation');
hold off;




% Calcul du nombre total d'evenement pour les differents cas

%calcul du nombre total d'evenements a partir de gamobs et l'exposure
nombre1=exposure*gamobs;
%calcul du nombre total d'evenements a partir de la simulation pour une largeur de 1 jour
nbbin1=temax;
for j =1:nbbin1;
bords1(j)=temax/nbbin1*(j-1);
end;
i=find(teobs<temax); %c'est le meme tableau pour les trois (quatres) largeurs
hn1=histc(teobs(i),bords1)*gamobs*exposure/length(teobs(i));
diff=length(teobs)-length(i); % c'est la meme valeur pour tous
nombre2=sum(hn1);
%calcul du nombre total d'evenements a partir de la simulation pour une largeur de 2 jours
nbbin2=temax/2;     
for j =1:nbbin2;
bords2(j)=temax/nbbin2*(j-1);
end;
hn2=histc(teobs(i),bords2)*gamobs*exposure/length(teobs(i));
nombre3=sum(hn2);
%calcul du nombre total d'evenements a partir de la simulation pour une largeur de 5 jours     
nbbin5=temax/5;
for j =1:nbbin5;
bords5(j)=temax/nbbin5*(j-1);
end;
hn5=histc(teobs(i),bords5)*gamobs*exposure/length(teobs(i));
nombre4=sum(hn5);
%calcul du nombre total d'evenements a partir de la simulation pour une largeur de 0.5 jour
nbbin=2*temax;
for j =1:nbbin;
bords(j)=temax/nbbin*(j-1);
end;
hn=histc(teobs(i),bords)*gamobs*exposure/length(teobs(i));
nombre5=sum(hn);
%calcul du nombre total d'evenements a partir des donnes experimentales
nombre6=length(teff);

%Affichage de ces resultats
disp(['Il n''y a pas eu prise en compte de ' num2str(diff) 'elements']);
disp(['nb d''evt avec gamobs :' num2str(nombre1) ]);
disp(['nb d''evt avec 1 j :' num2str(nombre2) ]);
disp(['nb d''evt avec 2 j :' num2str(nombre3) ]);
disp(['nb d''evt avec 5 j :' num2str(nombre4) ]);
disp(['nb d''evt avec 0.5 j :' num2str(nombre5) ]);
disp(['nb d''evt par l''experience :' num2str(nombre6) ]);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Trace dgamma (c'est le dgamma global et non pas celui apres le Monte-Carlo)
% Pour avoir quelque chose de coherent, trace plutot dgammaccepte
% figure(22);
% mgam=max(dgamma);
% y=0:mgam/100000:mgam;
% h4=histc(dgamma,y);
% h4=h4+1;
% h4=log(h4);
% bar(y,h4);
% title('distrib de dgamma');
% hold off;
%trace dgamma avec une echelle log pour la lisibilite
%figure(23);
%%dgamma=dgamma*1000+1;
%ldg=log(dgamma)/log(1000);
%hist(ldg,100000);
%title('distrib de dgamma echelle log1000');
%hold off;
%n'a pas trop d'interet




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Densite (Avec les espacements pour ramener les dimensions)

%figure(24);
%l=length(teff);
%milieu(1)=0;
%for i=2:l
%%    milieu(i)=(teff(i-1)+teff(i))/2;
%end;
%for i=1:(l-1)
%    ecart(i)=milieu(i+1)-milieu(i);
%end;
%ecart(l)=(teff(l)-milieu(l))*2;
%h3=1./ecart;
%bar(t,h3,'g');
%title('idee bizarre');
%hold off;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% En essayant de comparer des "resultats" normalises
% Ne sert a rien s'il n'y a pas une coherence avant
%figure(25);
%m1=max(h);
%h1=h1/m1;
%plot(centre,h1,'b-');
%hold on;
%h2=histc(teff,bords);
%m2=max(h2);
%h2=h2/m2;
%plot(centre,h2,'r-.');  ne pas tracer avec celui la
%bar(centre,h2,'g');
%title('comparaison en normalisee');
%hold off;




hold off;

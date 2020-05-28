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
% l = 4 *pi/180;    % direction d'observation en radian
% b = -1 *pi/180;

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

%----------------------------
% Param�tre pour le blending 
%---------------------------

f = 0.5; % fraction des évenements concernés par le blending
nbar = 1.257; % P(n) = fonction(nbar) = f avec P(n) la proba d'avoir n étoiles dans DeltaS


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



%------------------------
% Application du blending
%------------------------

% Preparation de la table pour le tirage aleatoire des luminosités


ll = (0:1e-5:1).*(Lsup-Linf)+Linf;
fll = probal(ll);
ifll = cumsum(fll);	% primitive 
ifll = ifll-ifll(1);
ifll = real(ifll./ifll(end));	% on fait en sorte que la primitive varie de 0 a 1

% Tirage est évenements concernés

ra = rand(size(te)); 
out = find(ra-f> 0);
in = find(ra-f<= 0);

% Tirage des luminosités (donc des flux)

ra1 = rand(size(te(in)));
ra2 = rand(size(te(in)));

flux1 = interp1(ifll,ll,ra1);
flux2 = interp1(ifll,ll,ra2);

B = flux1 ./ (flux1 + flux2);

% Blending

%opération pour pouvoir récupérer les infos plus tard : on conserve les indices et le résultats du blending par rapport à te

blend = zeros(size(in));

At=ampli(uT); %on utilse le seuil défini par l'expérience
Umin = rand(size(in));
Uobs= zeros(size(in));

%Calcul de Uobs
Uobs=maxampli(B.*ampli(Umin)+1-B);

%Calcul du facteur d'amplification
fact = sqrt((maxampli((At - 1)./B +1)).^2-Umin.^2)./(sqrt(1-Uobs.^2));
il=find(fact~=real(fact)); 
fact(il)=1; % donne parfois des nombres complexes si trop proche de l'amplification infinie, on pose donc une amplification = 1
	    % On simule te et on observe te,obs, donc nous voulons tracer te,obs

gmean = mean(fact);

%Récupération des résultats
teblend = te;
teblend(in)=te(in).*fact; % on a appliqué le blending à te et on a conservé l'ordre de te (important pour le blending après efficacité)


taurblend=taur * gmean * (nbar/(1-exp(-nbar)));
taurblend=real(taurblend);
disp(['tau avec blending (Alibert 2005)  = ' num2str(taurblend)]);

%-----------
%Choix de l'expérience à analyser
%------------

% exp_ogle_2006
% exp_macho_2005
exp_eros_2006

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
% affichage des resultats
%------------------------


%telechargement de la courbe du modèle

load ../graph/evenements_1.txt
te_model = evenements_1(:,5);

%Paramètre graph
nbre_bin = 100;
bin_max = 100;


%Trace la distribde te normalisée
[hist, edges] = histcounts(te, nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');
[hist_model, edges] = histcounts(te_model, nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');

%tracé de la distribution avec l'exposition 
[hist_ogle, edges] = histcounts(te, nbre_bin, 'BinLimits',[0,bin_max]);
[hist_ogle_model, edges] = histcounts(te_model, nbre_bin, 'BinLimits',[0,bin_max]);

%tracé distribution avec blending
[histb, edges] = histcounts(teblend, nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');

%Courbe expérimentale :
[hist_ogle_obs, edges] = histcounts(teobs, nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');
[hist_ogle_obs_b, edges] = histcounts(teobsblend, nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');


i=find(te<30);
i_model = find(te_model<30);
[hist_1, edges] = histcounts(te(i), nbre_bin, 'BinLimits',[0,bin_max], 'Normalization', 'probability');

centre = zeros(size(edges)-[0,1]);

for j =1:length(centre);
centre(j)=(edges(j)+edges(j+1))/2;
end

%Graph normalisé
figure(16)
hold on;
plot(centre, hist_1, 'black');
plot(centre, hist_model, 'red');
plot(centre, hist_ogle/length(te(i)), 'b')
title('Distribution de te normalisée');
xlabel('t_{e}')
ylabel('Nombre d''évènements par unité de t_{e}')

%graph en fonction de l'exposition
figure(17)
hold on;
plot(centre, hist_ogle_obs.*gamobs*exposure, 'red');
plot(centre, hist_ogle_obs_b*gamobs*exposure, 'black');
histogram(teff, nbre_bin, 'BinLimits',[0,bin_max])
title('Distribution de te vu par ogle');
xlabel('t_{e}')
ylabel('Nombre d''évènements par unité de t_{e}')


%graph noramlisé avec blending
passfigure(1);
hold on;
plot(centre, hist, 'red')
plot(centre, histb, 'black')
title('Blending black et sans blending rouge)');
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

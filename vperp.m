function vt = vperp(x,glr,glt,glp,sigrl,sigtl,sigpl,vrotl,L,gsr,gst,gsp,sigrs,sigts,sigps,vrots)

global cosb cosbl sinb L0 sinl vsp vsr vst vlr vlp vlt Ro


%------------------
% VITESSE DU SOLEIL
%------------------

vxsol=0; vzsol=0; vysol=-200e3;

%------------------
% VITESSE DU SOLEIL, Brunthaler et al. 2010
%------------------

% vxsol = -11.1; 
% vzsol = 7.25; 
% vysol = vrotdm(Ro,elev,0) + 12.24;

%-----------------------------------------------------------
% CALCUL DES COORDONNEES PAR RAPPORT AU CENTRE DE LA GALAXIE
%-----------------------------------------------------------

[R, z, th] = toGC(x.*L);
r  = sqrt(R.*R+z.*z);

%-------------------------------------------------------------------------
% Conversion des g en vitesse. g varie entre 0 et 1, et le v correspondant
% a une distribution gaussienne de centre 0 et de dispertion sig 
%-------------------------------------------------------------------------

vlr = sigrl.*erfinv(2.*glr-1);
vlt = sigtl.*erfinv(2.*glt-1);
vlp = sigpl.*erfinv(2.*glp-1)+vrotl;

%-----------------------------------------------
% calcul des angles pour conversion en cartesien
%-----------------------------------------------

sth = R./r;	cth = z./r;	sph = x.*L*cosb.*sinl./R;
cph = -(Ro-x.*L.*cosbl)./R;

%----------------------------------
% calcul de la vitesse en cartesien
%----------------------------------

vlx = vlr.*sth.*cph + vlt.*cth.*cph - vlp.*sph;
vly = vlr.*sth.*sph + vlt.*cth.*sph + vlp.*cph;
vlz = vlr.*cth - vlt.*sth;

%----------------------------------------------------------------------------
% on tient compte maintenant de la vitesse du Soleil et de la source (eq 3.21)
%-----------------------------------------------------------------------------

%----------------------------------------------
% Meme calcul mais pour la source cette fois-ci
%----------------------------------------------

[R, z, th] = toGC(L);	% utilisation des memes variables pour economiser
r  = sqrt(R.*R+z.*z);	% la memoire

vsr = sigrs.*erfinv(2.*gsr-1);
vst = sigts.*erfinv(2.*gst-1);
vsp = sigps.*erfinv(2.*gsp-1)+vrots;

sth = R./r;	cth = z./r;	sph = L*cosb.*sinl./R;
cph = -(Ro-L.*cosbl)./R;

vsx = vsr.*sth.*cph + vst.*cth.*cph - vsp.*sph;
vsy = vsr.*sth.*sph + vst.*cth.*sph + vsp.*cph;
vsz = vsr.*cth - vst.*sth;

%------------------------------
% CALCUL DE LA VITESSE RELATIVE
%------------------------------

vlx = vlx - (1-x).*vxsol + x.*vsx;
vly = vly - (1-x).*vysol + x.*vsy;
vlz = vlz - (1-x).*vzsol + x.*vsz;

%----------------------------------------------
% Vitesse projet�e le long de la ligne de vis�e
%----------------------------------------------

vr = cosbl.*vlx + cosb.*sinl.*vly +sinb.*vlz;
v  = vlx.*vlx + vly.*vly + vlz.*vlz;

%------------------------------------------------------------
% Norme de la vitesse perpendiculairement � la ligne de vis�e
%------------------------------------------------------------

vt = sqrt(v-vr.*vr);

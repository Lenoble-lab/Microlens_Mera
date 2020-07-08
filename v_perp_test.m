%------------------
% VITESSE DU SOLEIL
%------------------

vxsol=0; 
vysol=10;
vzsol=0; 

th = 0.01;
x = 0.9;

l = 0;
b = 0;

sinb = abs(sin(b));		cosb = cos(b);		cosl = cos(l);
cosbl=cos(b)*cos(l);		sinl = sin(l);
%-----------------------------------------------------------
% CALCUL DES COORDONNEES PAR RAPPORT AU CENTRE DE LA GALAXIE
%-----------------------------------------------------------

% [R, z, th] = toGC(x.*L);
%     r  = sqrt(R.*R+z.*z);

%-------------------------------------------------------------------------
% Conversion des g en vitesse. g varie entre 0 et 1, et le v correspondant
% a une distribution gaussienne de centre 0 et de dispertion sig 
%-------------------------------------------------------------------------

% vlr = sigrl.*erfinv(2.*glr-1).*0;
% vlt = sigtl.*erfinv(2.*glt-1).*0+vrotl;
% vlz = sigzl.*erfinv(2.*glz-1).*0;

vlr = 0;
vlt = 10;
vlz = 0;

%-----------------------------------------------
% calcul des angles pour conversion en cartesien
%-----------------------------------------------

%     sth = x.*L*cosb.*sinl./R;
%     cth = -(Ro-x.*L.*cosbl)./R;

cth = cos(th);
sth = sin(th);

%----------------------------------
% calcul de la vitesse en cartesien
%----------------------------------

vlx = vlr.*cth - vlt.*sth;
vly = vlr.*sth + vlt.*cth;

%-----------------------------------------------------------------------------
% on tient compte maintenant de la vitesse du Soleil et de la source (eq 3.21)
%-----------------------------------------------------------------------------

%----------------------------------------------
% Meme calcul mais pour la source cette fois-ci
%----------------------------------------------

% [R, z, th] = toGC(L);	% utilisation des memes variables pour economiser
%     r  = sqrt(R.*R+z.*z);	% la memoire

% vsr = sigrs.*erfinv(2.*gsr-1).*0;
% vst = sigts.*erfinv(2.*gst-1).*0+vrots;
% vsz = sigzs.*erfinv(2.*gsz-1);

vsr = 0;
vst = 10;
vsz = 0;

%     sth = x.*L*cosb.*sinl./R;
%     cth = -(Ro-x.*L.*cosbl)./R;

cth = cos(th);
sth = sin(th);
%----------------------------------
% calcul de la vitesse en cartesien
%----------------------------------

vsx = vsr.*cth - vst.*sth;
vsy = vsr.*sth + vst.*cth;


%------------------------------
% CALCUL DE LA VITESSE RELATIVE
%------------------------------

vlx = vlx - (1-x).*vxsol - x.*vsx;
vly = vly - (1-x).*vysol - x.*vsy;
vlz = vlz - (1-x).*vzsol - x.*vsz;

%----------------------------------------------
% Vitesse projet�e le long de la ligne de vis�e
%----------------------------------------------

vr = cosbl.*vlx + cosb.*sinl.*vly +sinb.*vlz;
v  = vlx.*vlx + vly.*vly + vlz.*vlz;

%------------------------------------------------------------
% Norme de la vitesse perpendiculairement � la ligne de vis�e
%------------------------------------------------------------

vt = sqrt(v-vr.*vr)

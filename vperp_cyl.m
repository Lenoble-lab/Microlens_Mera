function vt = vperp_cyl(x,glr,glt,glz,sigrl,sigtl,sigzl,vrotl,L,gsr,gst,gsz,sigrs,sigts,sigzs,vrots)

    global cosb cosbl sinb L0 sinl vsp vsr vst vlr vlp vlt Ro
    
    
    %------------------
    % VITESSE DU SOLEIL
    %------------------
    
    vxsol=0; vzsol=0; vysol=-200e3;
    
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
    vlt = sigtl.*erfinv(2.*glt-1)+vrotl;
    vlz = sigzl.*erfinv(2.*glz-1);
    
    %-----------------------------------------------
    % calcul des angles pour conversion en cartesien
    %-----------------------------------------------
    
%     sth = x.*L*cosb.*sinl./R;
%     cth = -(Ro-x.*L.*cosbl)./R;

    cth = -cos(th);
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
    
    [R, z, th] = toGC(L);	% utilisation des memes variables pour economiser
    r  = sqrt(R.*R+z.*z);	% la memoire
    
    vsr = sigrs.*erfinv(2.*gsr-1);
    vst = sigts.*erfinv(2.*gst-1)+vrots;
    vsz = sigzs.*erfinv(2.*gsz-1);
    
%     sth = x.*L*cosb.*sinl./R;
%     cth = -(Ro-x.*L.*cosbl)./R;

    cth = -cos(th);
    sth = sin(th);
    %----------------------------------
    % calcul de la vitesse en cartesien
    %----------------------------------
    
    vsx = vsr.*cth - vst.*sth;
    vsy = vsr.*sth + vst.*cth;
    
    
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
    
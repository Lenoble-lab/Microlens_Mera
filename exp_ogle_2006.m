%-----------------------------------
% Données OGLE II 2006 : Sumi et Al.
%-----------------------------------


%Exposure  OGLE II 2006 p 249
exposure = 1084267/10^(6)*1330/365.25;


ttobs=tau/(gam/1e6/365.25);
disp(['<tobs> (en jours) = ' num2str(ttobs)]);
N=gam*exposure;
disp(['nb d''evt  = ' num2str(N)]);

taur=gam*pi/2*uT*mean(te)/365.25/1e6;
taur=real(taur);
disp(['tau (avec gamma integré par MC) = ' num2str(taur)]);


%-----------------
%Evenements
%-----------------

teogle2006= [  44.7 ; 38.8; 46.3; 15.6; 41.2; 25.5; 7.5; 23.6; 22.9; 5.9; 24.2; 14.4; 16.0 ; 19.5; 65.5; 31.5; 18.9; 13.0; 53.8; 153.5; 140.9; 57.5; 66.8; 33.9; 17.2; 66.3; 12.0; 7.3; 33.9; 6.9; 14.2; 17.7];


effogle2006 = [  0.0209  ;  0.0376   ; 0.0528  ;  0.0888    ;0.1169  ;  0.1434  ;  0.2021 ;   0.2563   ; 0.3045  ;  0.3588   ; 0.3657   ; 0.4052   ; 0.4527   ; 0.4558;    0.4771;    0.5148  ;  0.5624  ;  0.6051   ; 0.6561 ;   0.6588  ;  0.6890    ;0.6905   ; 0.6892   ; 0.6903    ;0.6121    ;0.4710   ; 0.0433];
teffogle2006 = [1.1220;    1.4125;    1.7783;    2.2387;    2.8184;    3.5481;    4.4668;    5.6234;    7.0795;    8.9125;   11.2202;   14.1254;   17.7828;   22.3872;   28.1838;   35.4813;   44.6684;   56.2341;   70.7946;   89.1251;  112.2018;  141.2538;  177.8279;  223.8721;  281.8383;  354.8134;  446.6836];


%Données
teff = teffogle2006;
eff = effogle2006;

%-----------------------------------------------------------------------------------------------
% Interpolation lineaire de l'efficacite pour determiner la probabilite qu'un evt a d'etre garde
%-----------------------------------------------------------------------------------------------

% OGLE II

tinterpogle = teffogle2006;
effinterpogle = effogle2006;

teffmaxe=max(tinterpogle);
teffmine=min(tinterpogle);

i1 = find((te<=teffmaxe)&(te>=teffmine));
effsimogle = zeros(1,length(te));	% applique une efficacite nulle aux durees superieures et inferieures
effsimogleblend = zeros(1,length(te));
effsimogle(i1) = interp1(tinterpogle,effinterpogle,te(i1));
effsimogleblend(i1) = interp1(tinterpogle,effinterpogle,teblend(i1));


%-----------------------------------------------------------------------------
% Tirage d'un nombre aleatoire qui servira a decider si l'evt est garde ou non
%-----------------------------------------------------------------------------

raogle = rand(1,length(te))*max(effinterpogle);

%--------------------------------------------------------------------------------------------------------------------------
% compare le nombre aleatoire precedent a l'efficacite que l'on vient de calculer afin de decider si l'evt est garde ou non
%--------------------------------------------------------------------------------------------------------------------------

% On choisit l'efficacité ici en prenant les bons indices i


i = find(raogle-effsimogle<=0);

teobs = te(i);

ib = find(raogle-effsimogleblend<=0);

teobsblend = teblend(ib); % On récupère les éléments qui sont soumis au blending avec le calcul d'avant

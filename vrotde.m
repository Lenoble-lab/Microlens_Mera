% Courbe de rotation pour les etoiles du disque
%
% Parametre : les coordonnees galactiques
% retourne la vitesse vphi en m/s

function v=vlrot(R,z,th)
%v=ones(size(R)).*220e3;  	% rotation constante


% v=ones(size(R)).*200e3;   % utilise par E&G

%--------------------
%Brunthalter et al, 2010
%---------------------
% global Ro
% v_rot_sol = 239e3;
% 
% v = v_rot_sol *(1.00767 * (R./Ro).^0.0394 + 0.00712); 

v = vrotdm(R,z,th);
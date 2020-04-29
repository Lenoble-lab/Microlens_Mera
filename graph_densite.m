
clear
global dsup dinf Ro elev Rcoro

Rcoro = 3500;
elev = 0;
Ro = 8500;
dsup = 10000.;
dinf = 800.; %distance en parsec

global l b

% definition de la fenetre de Baade dans la majeure partie des articles :  l = 1 et b = -4
% A priori c'est cette definition qui est juste.
l = 1 *pi/180;    % direction d'observation en radian
b = -4 *pi/180;


global sinb cosb  cosbl sinl cosl
sinb = abs(sin(b));		cosb = cos(b);		cosl = cos(l);
cosbl=cos(b)*cos(l);		sinl = sin(l);


x = (0:1e-5:1).*(dsup-dinf)+dinf+Rcoro+5000;

[R, z, th] = toGC(x);

dens_zhao = rhozhao(R, z, th);
dens_dm = rhodm(R, z, th);
dens_de = rhode(R, z, th);

% figure(1)
% hold on;
% plot(x, dens_zhao);
% plot(x, dens_dm);
% plot(x, dens_de);
% % plot(x, R./dsup)
% % plot(x, z/350)
% 
% legend('bulbe (Zhao)', 'disque mince', 'disque épais');
% xlabel('distance au soleil (en pc))');
% ylabel('densité de masse en M_{sol}/pc^{3}');
% 
% format long
% m_bu = integral3(@rhozhao,-1000, 1000,-1000, 1000, 0, pi/4, 'method','iterated');
% disp(['masse bulbe ', num2str(m_bu*8)])

disp = sigtdm(R, z, th);

figure(2)
plot(R, disp);

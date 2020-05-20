
clear *
global dsup dinf Ro elev Rcoro

Rcoro = 3500;
elev = 0;
Ro = 8000;
dsup = 8500.;
dinf = 100.; %distance en parsec

global l b

% definition de la fenetre de Baade dans la majeure partie des articles :  l = 1 et b = -4
% A priori c'est cette definition qui est juste.
l = 1.5 *pi/180;    % direction d'observation en radian
b = -2.68 *pi/180;


global sinb cosb  cosbl sinl cosl
sinb = abs(sin(b));		cosb = cos(b);		cosl = cos(l);
cosbl=cos(b)*cos(l);		sinl = sin(l);


x = (0:1e-5:1).*(dsup-dinf)+dinf;

[R, z, th] = toGC(x);

% dens_zhao = rhozhao(R, z, th);
% dens_G2 = rhodwek(R, z, th);
% dens_E2 = rhostanek(R, z, th);
% dens_HetG = rhobuHetG(R, z, th);

% dens_d = rhodHetG(R, z, th);
% dens_dm = rhodm(R, z, th);
% dens_de = rhode(R, z, th);

i0 = find( R <= Rcoro );   
i1 = find( R > Rcoro);

vrot(i1) = vrotdm(R(i1),z(i1),th(i1));
vrot(i0) = vrotb(R(i0),z(i0),th(i0));

figure(1)
hold on;
plot(x,vrotdm(x,z,th).*1e-3);
title('Vitesse de rotation en fonction de R');
xlabel('distance au centre galactique')
ylabel('vitesse en km/s')
% plot(x, vrotb(R,z,th).*1e-2);
% plot(x, dens_zhao);
% plot(x, dens_E2);
% plot(x, dens_G2);
% legend('bulbe (Zhao)', 'E2 (Stanek)', 'G2 (dwek)', 'H&G');


% plot(x, dens_d)
% plot(x, dens_dm)
% plot(x, dens_de)  
% legend('H et G', 'dm', 'de')
% 
% xlabel('distance au soleil (en pc))');
% ylabel('densité de masse en M_{sol}/pc^{3}');

%% Integration pour le calcul de la masse
clear *
maxx = 1e5;
format long
m_bu = integral3(@rho_stanek_simple,0, maxx, 0, maxx, 0, maxx);
disp(['masse bulbe ', num2str(m_bu*8*1e-10)])

function rh = rho_stanek_simple(x,y,z)

x0=890;
y0=x0.*4.3/10;
z0=x0.*2.8/10;

M_b = 2e10;
rho0 = M_b/(x0*y0*z0*8*pi);

a = sqrt((x./x0).^2+(y./y0).^2+(z./z0).^2);

rh = rho0.*exp(-a);
end

function rh = rho_dweck_simple(X,Y,Z)
    
x0=890;
y0=x0.*4.3/10;
z0=x0.*2.8/10;
    
M_b = 1e10;
    rho0 = M_b/(6.57*pi*x0*y0*z0);
%     rho0 = 1/(2*pi*2^(3/2)*2*0.675*x0*y0*z0);
    sb2=sqrt(((X/x0).^2+(Y/y0).^2).^2+(Z/z0).^4);
    
    %-----------------
    % densite de masse
    %-----------------
    
    rh=rho0*(exp(-sb2/2)); 
end
function rh = rho_zhao_simple(X,Y,Z)
    
    x0=1490; 
    y0=580;
    z0=400;
    qa=0.6;
    
    M_b = 2;
    rho0 = M_b/0.9221;

    sb2=sqrt(((X/x0).^2+(Y/y0).^2).^2+(Z/z0).^4);
    sa=sqrt((qa*qa*(X.^2+Y.^2)+Z.^2)/(z0^2));

    %-------------  ----
    % densite de masse
    %-----------------

    rh=rho0*(exp(-sb2/2)+sa.^(-1.85).*exp(-sa));
end
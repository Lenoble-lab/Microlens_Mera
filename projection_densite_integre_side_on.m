%% test densité
clear; close all;
global X Z

X = 8500; Z = 0;
xx = (0:1e-2:1).*1e4;
zz = (0:1e-2:1).*3e3;

% [R, th] = meshgrid(R, th);
% % rho_bu = integral(@denssource(R, , th);
% % rho_disk = rhodHetG(R, zeros(size(R)), th);
% densite(R,th);
% [X,Y] = pol2cart(th, R);

dens_surf = zeros(size(xx,2), size(zz,2)); 

for k = 1:numel(xx)
    X = xx(k);
   
    if ( rem(k,10) == 0 )
disp(['compteur = ' num2str(k)]);

    end   
    for i = 1:numel(zz)
        Y = zz(i);
        dens_surf(k, i) = 2e6 * integral(@densite, 0, 1e5); %1e6 pour être en kpc
        
    end
end

[x_y, z_x] = meshgrid([-1 * fliplr(xx) xx],[-1 * fliplr(zz) zz]);
dens_surf = quadrant(dens_surf);


x_y = x_y * 1e-3; z_x = z_x * 1e-3;

figure(1)
% colormap(hot)
pcolor(x_y, z_x, log10(dens_surf)); 
shading interp;
hold on
colorbar;
caxis([6 10])
contour(x_y, z_x, log10(dens_surf), 'black', 'ShowText', 'on');


function res = quadrant(matrix)
%Construit une image symétrique à partir d'un unique quart d'image (éviter
%de faire trop de calcul)
    matrix = transpose(matrix);
    matrix = [fliplr(matrix) matrix];
    matrix = [flipud(matrix) ; matrix];    
    res = matrix;
end

function rh = densite(Y)
    
    global X Z Ro
    
    R = sqrt(X^2+Z^2);
    rh = zeros(size(R));
    %-------------
    %Bulbe
    %---------------
    
    x0=1580; %en parsec
    y0=620;
    z0=430;

    M_b= 1.8*1e10; 
    rho0 = M_b/(6.57*pi*x0*y0*z0);
    sb2=sqrt(((X/x0).^2+(Y/y0).^2).^2+(Z/z0).^4);
     
    rh = rh + rho0*(exp(-sb2/2)); 
    
    
    %-------------------
    % disque double expo
    %-------------------

    htd = 760;	        % Echelle de hauteur en pc
    Ltd = 3000; 	        % Echelle de longueur en pc
    rhtd = 0.05/20;        % densite au voisinage solaire, en Msol/pc^3

    rh = rh + rhtd.*exp(-R./Ltd-abs(Z)./htd+Ro./Ltd);

    %-------------------
    % disque double expo
    %-------------------

    hd=250;
    Ld=3000;
    rhsol=0.05;

    rh = rh + rhsol*exp(-R./Ld-abs(Z)./hd+Ro./Ld);
end


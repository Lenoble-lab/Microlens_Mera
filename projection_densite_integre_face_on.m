%% test densité
clear; close all;
global X Y
X = 8500;Y = 0;
xx = (0:1e-2:1).*1e4;
yy = (0:1e-2:1).*1e4;

%----------------------
%surface density selon les deux axes de la barre
%--------------------
dens_minor = zeros(size(yy));
dens_major = zeros(size(xx));

for i = 1:numel(yy)
    X = 0 ; Y = yy(i);
    dens_minor(i) = real(2e6 * integral(@densite, 0, 1e5));

end

for i = 1:numel(xx)
    X = xx(i) ; Y = 0;
    dens_major(i) = real(2e6 * integral(@densite, 0, 1e5));

end

xx_1 = [fliplr(xx)*-1 xx];
dens_minor = [fliplr(dens_minor) dens_minor];
dens_major = [fliplr(dens_major) dens_major];


figure(2)
semilogy(xx_1*1e-3, dens_minor)
hold on
semilogy(xx_1*1e-3, dens_major)
legend('Axe mineur', 'Axe majeur')
xlabel('X,Y (kpc)')
ylabel('\Sigma_{*} (M_{sun}.kpc^{-2})')


% [R, th] = meshgrid(R, th);
% % rho_bu = integral(@denssource(R, , th);
% % rho_disk = rhodHetG(R, zeros(size(R)), th);
% densite(R,th);
% [X,Y] = pol2cart(th, R);

dens_surf = zeros(size(xx,2), size(yy,2)); 

for k = 1:numel(xx)
    X = xx(k);
   
    if ( rem(k,10) == 0 )
disp(['compteur = ' num2str(k)]);

    end   
    for i = 1:numel(yy)
        Y = yy(i);
        dens_surf(k, i) = 2e6 * integral(@densite, 0, 1e5);
        
    end
end

[x_y, y_x] = meshgrid([-1 * fliplr(xx) xx],[-1 * fliplr(yy) yy]);
dens_surf = quadrant(dens_surf);


x_y = x_y * 1e-3; y_x = y_x * 1e-3;

figure(1)
% colormap(hot)
pcolor(x_y, y_x, log10(dens_surf)); 
shading interp;
hold on
colorbar;
caxis([7 10])
contour(x_y, y_x, log10(dens_surf), 'black', 'ShowText', 'on')
xlabel('X (kpc)')
ylabel('Y (kpc)')


function res = quadrant(matrix)
%Construit une image symétrique à partir d'un unique quart d'image (éviter
%de faire trop de calcul)
    matrix = transpose(matrix);
    matrix = [fliplr(matrix) matrix];
    matrix = [flipud(matrix) ; matrix];    
    res = matrix;
end



function rh = densite(z)
    %Densité face on (ne dépend que de z)
    global X Y Ro
    
    R = sqrt(X^2+Y^2);
    rh = zeros(size(R));
    %-------------
    %Bulbe
    %---------------
    
    x0=1580; %en parsec
    y0=620;
    z0=430;

    M_b= 1.8*1e10; 
    rho0 = M_b/(6.57*pi*x0*y0*z0);
    sb2=sqrt(((X/x0).^2+(Y/y0).^2).^2+(z/z0).^4);
     
    rh = rh + rho0*(exp(-sb2/2)); 
    
    
    %-------------------
    % disque double expo
    %-------------------

    htd = 760;	        % Echelle de hauteur en pc
    Ltd = 3000; 	        % Echelle de longueur en pc
    rhtd = 0.05/20;        % densite au voisinage solaire, en Msol/pc^3

    rh = rh + rhtd.*exp(-R./Ltd-abs(z)./htd+Ro./Ltd);

    %-------------------
    % disque double expo
    %-------------------

    hd=250;
    Ld=3000;
    rhsol=0.05;

    rh = rh + rhsol*exp(-R./Ld-abs(z)./hd+Ro./Ld);
end

% Densite de masse sous forme d'etoiles dans le disque
%
% Parametre : coordonnees cylindriques galactiques en pc R et z, et theta
% Retourne la densite en Msol/pc^3
%
%Formule de C. Han and A.P. Gould, Astrophys.J. 592 (2003) 172, astro-ph/0303309. (utilisé par Calchi Novatti)

function res = rhodHetG(R,z,th)

global Ro


H = 2.75e3 * Ro/8000;
h1 = 270 * Ro/8000;
h2 = 440 * Ro/8000;

eta = @(R) max([0.670.*ones(size(R)) ; 0.114 + R/(9.025e3)], [], 1);

%-----------------------%
%% disque simple uniquement %%
%-----------------------%
beta = 0;
rho_0 = 0.05 * eta(Ro);

%-----------------------%
%% disque épais et disque mince %%
%valeur de Han&Gould, réutilisée sans Iocco, 2018
%-----------------------%
% beta = 0.565;
% rho_0 = 4.93e-2;

res = rho_0./eta(R).*exp(-(R-Ro)./H).*((1-beta).*sech(z./(eta(R).*h1)).^2+beta.*exp(-abs(z)./(eta(R).*h2)));



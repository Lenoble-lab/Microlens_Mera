% Fonction a integrer pour le calcul de profondeur optique 

function res = dtau(x,L)
global Ro cosb cosbl sinb bet

res = rholens(x.*L).*x.*(1-x).*nsource(L).*L.^2;

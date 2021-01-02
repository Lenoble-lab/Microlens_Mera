%Fonction pour calculer gammax au début du MC
%Comme le calcul du max sur m et v est simple, on retire ces variables de
%la fonction et on optimise uniquement sur L et x
%On met en place la limitation sur dsup car sinon L tend vers +inf

function res = dgammmax(data)
x = data(1); L = data(2);

global dsup vlimit minf msup

if L<dsup
res = -rholens(x.*L).*sqrt(x.*(1-x)).*L.^1.5.*vlimit*sqrt(msup);	% Calcul du dgamma (cf eq 3.25 p. 57 de ma these)
else
    res = 0;
end
    %pas de terme nsource(L) car les distances
%sont tirees au hasard avec une loi proportionnelle
%a nsource



    
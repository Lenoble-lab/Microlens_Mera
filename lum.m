% relation magnitude/luminosité absolu : V - L : V = -2.5 log L + C
% grâce à la normalisation de la fl que nous allons effectuée, C n'a aucune importance et on peut prendre C = 0

function res = lum(V)

res = zeros(size(V));

res = 10.^(-V./(2.5));
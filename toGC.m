% [R, z, theta] = toGC(distance)
%
% Calcule les coordonnees cylindriques centrees sur le centre galactique
% d'un point situe a la distance specifiee du Soleil, dans une direction
% donnee. Ce qui revient a passer des coordonnees spheriques centrees
% sur le Soleil aux coordonnees cylindriques centrees sur la Galaxie.
% Les unites sont le pc partout.
%L'origine de theta et la direction du soleil

function [R, z, theta] = toGC(d)
global Ro sinl cosl cosbl cosb sinb elev
R = sqrt(Ro.*Ro+d.*d.*cosb.*cosb-2.*Ro.*d.*cosbl);
z = d.*sinb+elev;
theta = asin(d./R.*sinl.*cosb);
i = find(d*cosbl>=Ro); % cas des points situes au dela du centre galactique
if (length(i)>=1), theta(i) = pi-theta(i); end
    
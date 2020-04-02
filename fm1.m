% Fonction de masse dN/dm des lentilles.
% fonction FM1 de gilles

function pm = fm(m)
pm = zeros(size(m));

global minf msup

m1=1.0;


i1 = find(m>=minf & m<=m1);
if (length(i1)>=1)
  pm(i1) = 0.019*m(i1).^(-1.55);
end

i1 = find(m>=m1 & m<=msup);
if (length(i1)>=1)
  pm(i1) = 0.019*m(i1).^(-2.70);
end


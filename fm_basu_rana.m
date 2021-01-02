% Fonction de masse dN/dm des lentilles.
%IMF de basu et Rana (1992)

function pm = fm_basu_rana(m)
pm = zeros(size(m));

global minf msup

m1 = 0.08;
m2 = 0.5327;
m3 = 1.205;
m4 = 2;

i1 = find(m>=m1 & m<=m2);
i2 = find(m>=m2 & m<=m3);
i3 = find(m>=m3 & m<=m4);


if(length(i1)>=1)
  pm(i1) = 0.1919.*m(i1).^(-2.5);
end  

if(length(i2)>=1)
  pm(i2) = 0.1401 * m(i2).^(-3);
end  

if(length(i3)>=1)
  pm(i3) = 0.2861.*m(i3).^(-6.83);
end  
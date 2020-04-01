
% Fonction de masse dN/dm des lentilles.
% cette fonction de masse est celle de Chabrier PASP 2003 
% c'est l'equation 17 page 769

function pm = fmchab03(m)
pm = zeros(size(m));

global minf msup

m1=1.0;

sigma=0.69; % avec un ecart : 0.69 (-0.01/0.05)
m0=0.079;   % avec un ecart : 0.079 (-0.016/0.021)
c1=0.158/log(10);   % avec un ecart : 0.158 (0.051/-0.046)
c2=0.0441/log(10);  % la valeur 0.0441 vient de la continuite de la fonction de masse en 1 mase solaire
% c'est le resultat de 0.158*exp(-(log(0.079))^2/(2*069^2)) = 0.044096
alpha2=-1.35;       % attention c'est alpha = -1.35 +/- 0.3

i1 = find(m>=minf & m<=m1);
i2 = find(m>m1 & m<=msup);

if(length(i1)>=1)
  pm(i1) =c1*exp(-((log10(m(i1))-log10(m0)).^2)./(2*sigma^2))./m(i1);
end  

if(length(i2)>=1)
  pm(i2) =c2*((m(i2)).^(alpha2))./m(i2);
end  




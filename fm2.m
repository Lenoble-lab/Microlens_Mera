% Fonction de masse dN/dm des lentilles.

function pm = fm2(m)
pm = zeros(size(m));

global minf msup

sigma=0.627;
m0=0.1;

i1 = find(m>=minf & m<=msup);
if(length(i1)>=1)
  pm(i1) =0.141*log(10)*exp(-((log10(m(i1))-log10(m0)).^2)./(2*sigma^2))./m(i1);
end  




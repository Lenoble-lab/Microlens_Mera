% Fonction de masse dN/dm des lentilles.

function pm = fm3(m)
pm = zeros(size(m));

global minf msup


A=3.0;
mbar=716.4;
beta=0.25;
alpha=3.3;


i1 = find(m>=minf & m<=msup);
if (length(i1)>=1)
  pm(i1) = A.*exp(-(mbar./m(i1)).^beta).*(m(i1).^(-alpha));
end



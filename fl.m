% fonction de luminosité dN/dL cf Alibert 2005 + Holtzman 1998 (cité dans Alibert) dN/dL proportionnelle à L^-2

function res = fl(L)

res = zeros(size(L));

global Vinf Vsup

Vinf = 22;
Vsup = 16;

Linf=lum(Vinf);
Lsup=lum(Vsup);

il=find((L<=Lsup)&(L>=Linf));

res(il) = L(il).^(-2);

res;

%IMF de l'article de chabrier 2014
%Sur la fonction de masse au centre de la galaxie 

function pm = fmchab_2014_cas_1(m)
pm = zeros(size(m));

global minf msup

%Cas 1
x = 1.35;
m0 = 0.35;
n_c = 14;
sigma = 0.607;
m_c = 0.025;
A_h = 0.417;

i1 = find(m>=minf & m<=m0);
i2 = find(m>m0 & m<=msup);
A_l = A_h*n_c^(x/2);

if(length(i1)>=1)
  pm(i1) =A_l*m0^(-x) * exp(-((log10(m(i1))-log10(m_c)).^2)./(2*sigma^2)).*m(i1).^(-1);
end  

if(length(i2)>=1)
  pm(i2) =A_h*((m(i2)).^(-x))./m(i2);
end
end
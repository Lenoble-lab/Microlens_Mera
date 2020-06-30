function pm = fmchab05(m)
pm = zeros(size(m));

global minf msup

m1=1.0;

sigma=0.55;
m0=0.2;
c1=0.093/log(10);
c2=0.0415/log(10);
alpha2=-1.35;       % attention c'est alpha = -1.35 +/- 0.3

i1 = find(m>=minf & m<=m1);
i2 = find(m>m1 & m<=msup);

if(length(i1)>=1)
  pm(i1) =c1*exp(-((log10(m(i1))-log10(m0)).^2)./(2*sigma^2))./m(i1);
end  

if(length(i2)>=1)
  pm(i2) =c2*((m(i2)).^(alpha2))./m(i2);
 
end  
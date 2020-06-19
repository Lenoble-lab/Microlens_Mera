%fm de Kroupa de 2001

function pm = fm(m)
pm = zeros(size(m));
    
global minf msup

m1 = 0.01;
m2 = 0.08;
m3 = 0.5

i1 = find(m>=m1 & m<=m2 & m>=minf);
i2 = find(m>=m2 & m<=m3);
i3 = find(m>=m3 et m<msup);

c = 1

if (length(i1)>=1)
    pm(i1) = c * m(i1).^(-0.3);
end

if (length(i2)>=1)
    pm(i2) = c * 0.08^(-0.3)/0.08^(-1.3) * m(i2).^(-1.3);
end

if (length(i1)>=1)
    pm(i3) = c * 0.08^(-0.3)/0.08^(-1.3) * 0.5^(-1.3)/0.05^(-2.3)*m(i3).^(-2.3);
end
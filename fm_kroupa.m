%fm de Kroupa de 2001

function pm = fm_kroupa(m)
pm = zeros(size(m));
    
global minf msup

m1 = 0.01;
m2 = 0.08;
m3 = 0.5;

i1 = find(m>=m1 & m<m2 & m>=minf);
i2 = find(m>=m2 & m<m3);
i3 = find(m>=m3 & m<msup);

c = 1/0.1405;

alpha_bd = -0.3;
alpha_MS = -1.3;
alpha = -2.3;


if (length(i1)>=1)
    pm(i1) = c * m(i1).^(alpha_bd);
end

if (length(i2)>=1)
    pm(i2) = c * 0.08^(alpha_bd)/0.08^(alpha_MS) * m(i2).^(alpha_MS);
end

if (length(i1)>=1)
    pm(i3) = c * 0.08^(alpha_bd)/0.08^(alpha_MS) * 0.5^(alpha_MS)/0.5^(alpha)*m(i3).^(alpha);
end
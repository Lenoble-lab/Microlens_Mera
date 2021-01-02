%IMF du modèle de Besancon 
%Awiphan, 2016
%Fonction du disque

function pm = fm_awiphan_disk(m)
pm = zeros(size(m));

global minf msup

% Thin disk
m1 = 0.01;
m2 = 0.079;
m3 = 1;



c = 1;


%--------------------
%Thin disk
%------------------
alpha_bd = -0.4;
alpha_MS = -1.6;
alpha = -3;


i1 = find(m>=m1 & m<m2 & m>=minf);
i2 = find(m>=m2 & m<m3);
i3 = find(m>=m3 & m<msup);

if (length(i1)>=1)
    pm(i1) = c * m(i1).^(alpha_bd);
end

if (length(i2)>=1)
    pm(i2) = c * m2^(alpha_bd)/m2^(alpha_MS) * m(i2).^(alpha_MS);
end

if (length(i3)>=1)
    pm(i3) = m2^(alpha_bd)/m2^(alpha_MS) * m3^(alpha_MS)/m3^(alpha)*m(i3).^(alpha);
end 
%fm de Kroupa de 2001

function pm = fm_kroupa(m)
pm = zeros(size(m));
    
global minf msup

m1 = 0.01;
m2 = 0.08;
m3 = 0.5;



c = 1/0.1405;


%--------------------
%Kroupa 2001
%------------------
alpha_bd = -0.3;
alpha_MS = -1.3;
alpha = -2.3;

%-----------------------
%Valeurs de Calchi-Novatti
%-------------------------
% alpha_bd = -1.6;
% alpha_MS = -1.7;
% alpha = -2;

% alpha_bd = 1;
% alpha_MS = -0.15;
% alpha = -1;




i1 = find(m>=m1 & m<m2 & m>=minf);
i2 = find(m>=m2 & m<m3);
i3 = find(m>=m3 & m<msup);

if (length(i1)>=1)
    pm(i1) = c * m(i1).^(alpha_bd);
end

if (length(i2)>=1)
    pm(i2) = c * m2^(alpha_bd)/m2^(alpha_MS) * m(i2).^(alpha_MS);
end

if (length(i1)>=1)
    pm(i3) = c * m2^(alpha_bd)/m2^(alpha_MS) * m3^(alpha_MS)/m3^(alpha)*m(i3).^(alpha);
end
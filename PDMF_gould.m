
%Updating for a 10Gyr pop according to Gould 2000

function pm = PDMF_gould(pm,m)
i_WD = find(m>1 & m<8);
m_WD = 0.6;
[~,idx] = min(abs(m-m_WD));
pm(idx) = pm(idx) + sum(pm(i_WD));
% disp(sum(pm(i_WD)));

pm(i_WD) = zeros(size(i_WD));

i_NS = find(m>=8 & m<40);
m_NS = 1.35;
sig_NS = 0.04;
% [~,idx] = min(abs(m-m_NS));
% pm(idx) = pm(idx) + sum(pm(i_NS));
idx = find(abs(m-m_NS)<0.3);
pm(idx) = pm(idx) + sum(pm(i_NS))/(sig_NS * sqrt(2*pi)) * exp(-(m(idx)-m_NS).^2/(2*sig_NS^2));
pm(i_NS) = zeros(size(i_NS));

i_BH = find(m>40);
m_BH = 5;
sig_BH = 1;
% [~,idx] = min(abs(m-m_BH));
idx = find(abs(m-m_BH)<3);
pm(idx) = pm(idx) + sum(pm(i_BH))/(sig_BH * sqrt(2*pi)) * exp(-(m(idx)-m_BH).^2/(2*sig_BH^2));
% disp(sum(pm(i_BH)));

pm(i_BH) = zeros(size(i_BH)); 
end

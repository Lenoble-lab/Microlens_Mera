
%Updating for a 10Gyr pop according to Gould 2000

function pdmf = PDMF_gould(imf,m)

pm = zeros(size(m));

i_WD = find(m>1 & m<8);
m_WD = 0.6;
sig_WD = 0.10;
% [~,idx] = min(abs(m-m_WD));
% pm(idx) = pm(idx) + sum(pm(i_WD));
idx = find(abs(m-m_WD)<0.4);
pm(idx) = pm(idx) + sum(imf(i_WD))/(sig_WD * sqrt(2*pi)) * exp(-(m(idx)-m_WD).^2/(2*sig_WD^2));
disp(sum(pm(idx))/sum(imf(i_WD)));

i_NS = find(m>=8 & m<40);
m_NS = 1.35;
sig_NS = 0.04;
% [~,idx] = min(abs(m-m_NS));
% pm(idx) = pm(idx) + sum(pm(i_NS));
idx = find(abs(m-m_NS)<0.3);
pm(idx) = pm(idx) + sum(imf(i_NS))/(sig_NS * sqrt(2*pi)) * exp(-(m(idx)-m_NS).^2/(2*sig_NS^2));
disp(sum(pm(idx))/sum(imf(i_NS)));

i_BH = find(m>40);
m_BH = 5;
sig_BH = 1;
% [~,idx] = min(abs(m-m_BH));
idx = find(abs(m-m_BH)<3.3);
pm(idx) = pm(idx) + sum(imf(i_BH))/(sig_BH * sqrt(2*pi)) * exp(-(m(idx)-m_BH).^2/(2*sig_BH^2));
disp(sum(pm(idx))/sum(imf(i_BH)));

imf([i_WD i_NS i_BH]) = zeros(size([i_WD i_NS i_BH]));
pdmf = pm + imf;

end

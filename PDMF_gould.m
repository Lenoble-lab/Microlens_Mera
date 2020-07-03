
%Updating for a 10Gyr pop according to Gould 2000

function pdmf = PDMF_gould(imf,m)

pm = zeros(size(m));
pas = m(2)-m(1);

%WD
i_WD = find(m>1 & m<8);     %sélection des masses qui évolueront en WD
m_WD = 0.6;                 % Masse moy finale
sig_WD = 0.16;              %écart typer
% [~,idx] = min(abs(m-m_WD));
idx = find(abs(m-m_WD)<0.3);       %Selection de la planche sur laquelle on va projeter les données
if length(idx)>1
gauss_WD = @(x) 1/(sig_WD * sqrt(2*pi)) * exp(-(x-m_WD).^2/(2*sig_WD^2));   %Fonctino gaussienne (pour la normalisation)
pm(idx) = sum(imf(i_WD)) * pas * 1/(sig_WD * sqrt(2*pi)) * exp(-(m(idx)-m_WD).^2/(2*sig_WD^2))/integral(gauss_WD, m(idx(1)), max(m(idx)));
end
%Neutron star
i_NS = find(m>=8 & m<40);
m_NS = 1.35;
sig_NS = 0.04;
gauss_NS = @(x) 1/(sig_NS * sqrt(2*pi)) * exp(-(x-m_NS).^2/(2*sig_NS^2));

idx = find(abs(m-m_NS)<0.4);
if length(idx)>1
pm(idx) = pm(idx) + sum(imf(i_NS))*pas * gauss_NS(m(idx))/integral(gauss_NS, m(idx(1)), max(m(idx)));
end
%Black Hole
i_BH = find(m>40);
m_BH = 5;
sig_BH = 1;
gauss_BH = @(x) 1/(sig_BH * sqrt(2*pi)) * exp(-(x-m_BH).^2/(2*sig_BH^2));

% [~,idx] = min(abs(m-m_BH));
idx = find(abs(m-m_BH)<3.3);
if length(idx)>1
pm(idx) = pm(idx) + sum(imf(i_BH))*pas * 1/(sig_BH * sqrt(2*pi)) * exp(-(m(idx)-m_BH).^2/(2*sig_BH^2))/integral(gauss_BH, m(idx(1)), max(m(idx)));
end

disp([sum(imf(i_WD)), sum(imf(i_NS)), sum(imf(i_BH))]/sum(imf(m<1)))

imf([i_WD i_NS i_BH]) = zeros(size([i_WD i_NS i_BH]));
pdmf = pm + imf;

end


%Updating for a 10Gyr pop according to Gould 2000
%m : tirage aléatoire des masses selon l'imf considérée

function [m, frac_N, frac_M, frac_eve, star_pop] = PDMF_gould_1(m)

%BD
i_BD = find(m<0.07);
fract_BD = length(i_BD);
fract_M_BD = sum(m(i_BD));
fract_eve_BD = sum(sqrt(m(i_BD)));
star_pop(i_BD) = ones(size(i_BD));

%MS
i_MS = find(m<1 & m>0.07);
fract_MS = length(i_MS);
fract_M_MS = sum(m(i_MS));
fract_eve_MS = sum(sqrt(m(i_MS)));
star_pop(i_MS) = ones(size(i_MS)).*2;


%WD
i_WD = find(m>0.7 & m<8);     %sélection des masses qui évolueront en WD
m_WD = 0.6;                 % Masse moy finale
sig_WD = 0.04;              %écart type
% gauss_WD = @(x) 1/(sig_WD * sqrt(2*pi)) * exp(-(x-m_WD).^2/(2*sig_WD^2));   %Fonctino gaussienne (pour la normalisation)
fract_WD = length(i_WD);
m(i_WD) = sig_WD*randn(size(i_WD)) + m_WD;
fract_m_WD = sum(m(i_WD));
star_pop(i_WD) = ones(size(i_WD)).*3;


%Neutron star
i_NS = find(m>=8 & m<40);
m_NS = 1.35;
sig_NS = 0.04;
% gauss_NS = @(x) 1/(sig_NS * sqrt(2*pi)) * exp(-(x-m_NS).^2/(2*sig_NS^2));
fract_NS = length(i_NS);
m(i_NS) = sig_NS*randn(size(i_NS)) + m_NS;
fract_m_NS = sum(m(i_NS));
star_pop(i_NS) = ones(size(i_NS)).*4;


%Black Hole
i_BH = find(m>40);
m_BH = 5;
sig_BH = 1;
% gauss_BH = @(x) 1/(sig_BH * sqrt(2*pi)) * exp(-(x-m_BH).^2/(2*sig_BH^2));
fract_BH = length(i_BH);
m(i_BH) = sig_BH*randn(size(i_BH)) + m_BH;
fract_m_BH = sum(m(i_BH));
star_pop(i_BH) = ones(size(i_BH)).*5;

fract_WD = length(i_WD);
fract_m_WD = sum(m(i_WD));
fract_eve_WD = sum(sqrt(m(i_WD)));

fract_NS = length(i_NS);
fract_m_NS = sum(m(i_NS));
fract_eve_NS = sum(sqrt(m(i_NS)));

fract_BH = length(i_BH);
fract_m_BH = sum(m(i_BH));
fract_eve_BH = sum(sqrt(m(i_BH)));


frac_N = [fract_BD fract_MS fract_WD fract_NS fract_BH]./length(m);
frac_M = [fract_M_BD fract_M_MS fract_m_WD fract_m_NS fract_m_BH]/sum(m);
frac_eve = [fract_eve_BD fract_eve_MS fract_eve_WD fract_eve_NS fract_eve_BH]/sum([fract_eve_BD fract_eve_MS fract_eve_WD fract_eve_NS fract_eve_BH]);

end

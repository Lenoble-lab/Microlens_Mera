
%Updating for a 10Gyr pop according to Maraston 1998
%m : tirage aléatoire des masses selon l'imf considérée
function [m, frac_N, frac_M]  = PDMF_Maraston_1(m)

pm = zeros(size(m));

i_WD = find(m>1 & m<8.5);
m(i_WD) = 0.077.*m(i_WD)+0.48;


i_NS = find(m>=8.5 & m<40);
m(i_NS) = 1.4;


i_BH = find(m>40);
m(i_BH) = 0.5*m(i_BH);

fract_WD = length(i_WD);
fract_m_WD = sum(m(i_WD));

fract_NS = length(i_NS);
fract_m_NS = sum(m(i_NS));

fract_BH = length(i_BH);
fract_m_BH = sum(m(i_BH));

frac_N = [fract_WD fract_NS fract_BH]./length(m);
frac_M = [fract_m_WD fract_m_NS fract_m_BH]/sum(m);
end
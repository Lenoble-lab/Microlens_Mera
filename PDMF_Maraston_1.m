
%Updating for a 10Gyr pop according to Maraston 1998
%m : tirage aléatoire des masses selon l'imf considérée
%On stocke trois liste pour avoir les fractions en nombre, en masse et avec
%le facteur 1/sqrt(m) pour chaque type d'objet 
function [m, frac_N, frac_M, frac_eve, star_pop]  = PDMF_Maraston_1(m)

star_pop = zeros(size(m));

i_BD = find(m<0.07);
fract_BD = length(i_BD);
fract_M_BD = sum(m(i_BD));
fract_eve_BD = sum(sqrt(m(i_BD)));
star_pop(i_BD) = ones(size(i_BD));

i_MS = find(m<1 & m>0.07);
fract_MS = length(i_MS);
fract_M_MS = sum(m(i_MS));
fract_eve_MS = sum(sqrt(m(i_MS)));
star_pop(i_MS) = ones(size(i_MS)).*2;

i_WD = find(m>1 & m<8.5);
m(i_WD) = 0.077.*m(i_WD)+0.48;
star_pop(i_WD) = ones(size(i_WD)).*3;

i_NS = find(m>=8.5 & m<40);
m(i_NS) = 1.4;
star_pop(i_NS) = ones(size(i_NS)).*4;

i_BH = find(m>40);
m(i_BH) = 0.5*m(i_BH);
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

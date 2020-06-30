
%Updating for a 10Gyr pop according to Maraston 1998

function pdmf = PDMF_Maraston(imf,m)

pm = zeros(size(m));

i_WD = find(m>1 & m<8.5);
m_WD = 0.077.*m(i_WD)+0.48;
for m_WD_i =m_WD
[~,idx] = min(abs(m-m_WD_i));
pm(idx) = pm(idx) + imf(idx);
% disp(sum(pm(i_WD)));
end

i_NS = find(m>=8.5 & m<40);
m_NS = 1.4;
[~,idx] = min(abs(m-m_NS));
pm(idx) = pm(idx) + sum(imf(i_NS));


i_BH = find(m>40);
m_BH = 0.5*m(i_BH);
for m_BH_i =m_BH
[~,idx] = min(abs(m-m_BH_i));
pm(idx) = pm(idx) + imf(idx);
% disp(sum(pm(i_WD)));
end
imf([i_WD i_NS i_BH]) = zeros(size([i_WD i_NS i_BH]));
pdmf = pm + imf;
end

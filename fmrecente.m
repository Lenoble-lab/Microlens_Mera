% Fonction de masse dN/dm des lentilles.
% cette fonction de masse est la plus recente : 2005

% global minf msup
% 
% minf=0.01;
% msup=100;
% 
% m = (0:1e-5:1).*(msup-minf)+minf;
% 
% fm05 = fmrecente(m);
% figure(1)
% loglog(m, fm05)


function pm = fmrecente(m)
pm = zeros(size(m));

global minf msup

m1=1.0;

sigma=0.55;
m0=0.2;
c1=0.093/log(10);
c2=0.0415/log(10);
alpha2=-1.35;       % attention c'est alpha = -1.35 +/- 0.3

i1 = find(m>=minf & m<=m1);
i2 = find(m>m1 & m<=msup);
if(length(i1)>=1)
  pm(i1) =c1*exp(-((log10(m(i1))-log10(m0)).^2)./(2*sigma^2))./m(i1);
end  

if(length(i2)>=1)
  pm(i2) =c2*((m(i2)).^(alpha2))./m(i2);
end
% pm = fm_kroupa(m);
% pm = fmchab_modi(m);

% %Calculate de PDMF according to Gould 2000
% 
% i_WD = find(m>1 & m<8);
% m_WD = 0.6;
% [~,idx] = min(abs(m-m_WD));
% pm(idx) = pm(idx) + sum(pm(i_WD));
% % disp(sum(pm(i_WD)));
% 
% pm(i_WD) = zeros(size(i_WD));
% 
% i_NS = find(m>=8 & m<40);
% m_NS = 1.35;
% sig_NS = 0.04;
% % [~,idx] = min(abs(m-m_NS));
% % pm(idx) = pm(idx) + sum(pm(i_NS));
% idx = find(abs(m-m_NS)<0.3);
% pm(idx) = pm(idx) + sum(pm(i_NS))/(sig_NS * sqrt(2*pi)) * exp(-(m(idx)-m_NS).^2/(2*sig_NS^2));
% pm(i_NS) = zeros(size(i_NS));
% 
% i_BH = find(m>40);
% m_BH = 5;
% sig_BH = 1;
% % [~,idx] = min(abs(m-m_BH));
% idx = find(abs(m-m_BH)<3);
% pm(idx) = pm(idx) + sum(pm(i_BH))/(sig_BH * sqrt(2*pi)) * exp(-(m(idx)-m_BH).^2/(2*sig_BH^2));
% % disp(sum(pm(i_BH)));
% 
% pm(i_BH) = zeros(size(i_BH));
% 
% 

end

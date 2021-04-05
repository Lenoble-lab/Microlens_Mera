global minf msup

minf=0.01;
msup=100;

integral(@rapport_masse, 0.01, 0.07)/integral(@rapport_masse, 0.01, 100)
integral(@rapport_masse, 5, 50)/integral(@rapport_masse, 0.01, 100)

integral(@rapport_masse, 50, 100)

m = (0:1e-5:1).*(msup-minf)+minf;
fm05 = fmchab05(m)./integral(@fmchab05, minf, msup);
fm03 = fmchab03(m)./integral(@fmchab03, minf, msup);
%fm_maraston = PDMF_Maraston_1(m);
% fm_eff = rapport_masse(m)./integral(@rapport_masse, 0.01, 100);
%fm_gould = PDMF_gould_1(m);
% fm_gould = fm_basu_rana(m);
% fm_recente = fmrecente(m)/integral(@fmrecente, 0.01, 100);
fm_krou = fm_kroupa(m)./integral(@fm_kroupa, minf, msup);
fm_awi_disk = fm_awiphan_disk(m)./integral(@fm_awiphan_disk, minf, msup);
fm_awi_bulge = fm_awiphan_bulge(m)./integral(@fm_awiphan_bulge, minf, msup);
fm_cas_1 = fmchab_2014_cas_1(m)./integral(@fmchab_2014_cas_1, minf, msup);
fm_cas_2 = fmchab_2014_cas_2(m)./integral(@fmchab_2014_cas_2, minf, msup);
fm_MW = fmchab_2014_MW(m)./integral(@fmchab_2014_MW, minf, msup);

% fm_awi = fm_awiphan(m)./sum(fm_awiphan(m))/((msup-minf)/length(m));

[~,idx] = min(abs(m-1.4));

if ishandle(1)
    close(1)
end
figure(1);
set(gca, 'YScale', 'log')
hold on;
plot(log10(m),fm_cas_1);
plot(log10(m),fm_cas_2);
plot(log10(m),fm_MW);
plot(log10(m),fm03);
plot(log10(m),fm05);
% plot(log10(m),fm_krou);
% plot(log10(m),fm_awi_bulge);
% plot(log10(m),fm_awi_disk);
% plot(log10(m),fm_gould);
% plot(log10(m(idx)),fm_maraston(idx), 'x', 'Color', [0.8500 0.3250 0.0980]);
ylabel('\xi(log(m))')
xlabel('log_{10} (M_{sol})');
% legend('Chabrier 2005', 'chabrier 2003', 'Kroupa 2001', 'Awiphan, bulge', 'Awiphan, disk')
legend('cas 1', 'cas 2', 'cas MW', 'chabrier 2003', 'Chabrier 2005')
% axis([-2 2 1e-6 10]);



% fm05 = fmrecente(m)/integral(@fmrecente, 0.01, 100);
% fm03 = fmchab03(m)/integral(@fmchab03, 0.01, 100);
% fmkroupa = fm_kroupa(m)/integral(@fm_kroupa, 0.01, 100);
% fm_modi = fmchab_modi(m)/integral(@fmchab_modi, 0.01, 100);
% 
% figure(2);
% set(gca, 'YScale', 'log')
% hold on;
% plot(log10(m),fmkroupa);
% plot(log10(m),fm03);
% plot(log10(m),fm05, 'x');
% plot(log10(m),fm_modi);
% hold off
% 
% legend('MF1 (Kroupa)','MF2 (2003)','MF3 (2005)', 'fm chabier modif');
% ylabel('\xi(log(m))');
% xlabel('log_{10} (M_{sol})');
% %axis([0,100,0,0.05]);

% fm_bu = fmbu(m)/integral(@fmbu, 0.01, 100);
% fm_de = fm_bu;
% fm_dm = fm_bu;
% fm_de = fmde(m)/integral(@fmde, 0.01, 100);
% fm_dm = fmdm(m)/integral(@fmdm, 0.01, 100);
% 
% figure(2);
% set(gca, 'YScale', 'log')
% hold on;
% plot(log10(m),fm_bu);
% plot(log10(m),fm_de);
% plot(log10(m),fm_dm);
% hold off
% 
% legend('FM bulbe', 'FM de','FM dm');
% ylabel('\xi(log(m))');
% xlabel('log_{10} (M_{sol})');
%axis([0,100,0,0.05]);


function pm = fmchab_2014_MW(m)
pm = zeros(size(m));

global minf msup

% Cas MW
x = 1.35;
m0 = 2.0;
n_c = 11;
sigma = 0.589;
m_c = 0.18;
A_h = 0.649;

i1 = find(m>=minf & m<=m0);
i2 = find(m>m0 & m<=msup);
A_l = A_h*n_c^(x/2);

if(length(i1)>=1)
  pm(i1) =A_l*m0^(-x) * exp(-((log10(m(i1))-log10(m_c)).^2)./(2*sigma^2)).*m(i1).^(-1);
end  

if(length(i2)>=1)
  pm(i2) =A_h*((m(i2)).^(-x))./m(i2);
end
end
function pm = fmchab_2014_cas_2(m)
pm = zeros(size(m));

global minf msup

%cas 2
x = 1.6;
m0 = 0.35;
n_c = 11;
sigma = 0.531;
m_c = 0.032;
A_h = 0.390;

i1 = find(m>=minf & m<=m0);
i2 = find(m>m0 & m<=msup);
A_l = A_h*n_c^(x/2);

if(length(i1)>=1)
  pm(i1) =A_l*m0^(-x) * exp(-((log10(m(i1))-log10(m_c)).^2)./(2*sigma^2)).*m(i1).^(-1);
end  

if(length(i2)>=1)
  pm(i2) =A_h*((m(i2)).^(-x))./m(i2);
end
end
function pm = fmchab_2014_cas_1(m)
pm = zeros(size(m));

global minf msup

%Cas 1
x = 1.35;
m0 = 0.35;
n_c = 14;
sigma = 0.607;
m_c = 0.025;
A_h = 0.417;

i1 = find(m>=minf & m<=m0);
i2 = find(m>m0 & m<=msup);
A_l = A_h*n_c^(x/2);

if(length(i1)>=1)
  pm(i1) =A_l*m0^(-x) * exp(-((log10(m(i1))-log10(m_c)).^2)./(2*sigma^2)).*m(i1).^(-1);
end  

if(length(i2)>=1)
  pm(i2) =A_h*((m(i2)).^(-x))./m(i2);
end
end
function pm = rapport_masse(m)

pm = fmchab05(m);
% pm = fm_kroupa(m);
% pm = fmchab_modi(m);
pm = PDMF_Maraston(pm,m).*sqrt(m);
% pm = PDMF_gould(pm,m).*m;
end

function pm = fmrecente_maraston(m)

pm = fmchab05(m);
% pm = fm_kroupa(m);
% pm = fmchab_modi(m);
pm = PDMF_Maraston(pm,m);
% pm = PDMF_gould(pm,m);
end

function pm = fmrecente_gould(m)

pm = fmchab05(m);
% pm = fm_kroupa_modif(m);
% pm = fmchab_modi(m);
pm = PDMF_gould(pm,m);
end


global minf msup

minf=0.01;
msup=100;

m = (0:1e-5:1).*(msup-minf)+minf;
fm05 = fmchab05(m);
% fm_maraston = PDMF_Maraston(fm05,m);
fm_maraston = fmrecente_maraston(m);
fm_gould = fmrecente_gould(m);
% fm_recente = fmrecente(m)/integral(@fmrecente, 0.01, 100);

[~,idx] = min(abs(m-1.4));

if ishandle(1)
    close(1)
end
figure(1);
set(gca, 'YScale', 'log')
hold on;
plot(log10(m),fm05);
plot(log10(m),fm_maraston, 'm');
plot(log10(m),fm_gould);
plot(log10(m(idx)),fm_maraston(idx).*sqrt(m(idx)), 'xm');
ylabel('\xi(log(m))');
xlabel('log_{10} (M_{sol})');
legend('fm 005', 'fm Maraston', 'fm Gould')
axis([-1.5 1.5 1e-7 1]);



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

fm_bu = fmbu(m)/integral(@fmbu, 0.01, 100);
fm_de = fm_bu;
fm_dm = fm_bu;
fm_de = fmde(m)/integral(@fmde, 0.01, 100);
fm_dm = fmdm(m)/integral(@fmdm, 0.01, 100);

figure(2);
set(gca, 'YScale', 'log')
hold on;
plot(log10(m),fm_bu);
plot(log10(m),fm_de);
plot(log10(m),fm_dm);
hold off

legend('FM bulbe', 'FM de','FM dm');
ylabel('\xi(log(m))');
xlabel('log_{10} (M_{sol})');
% %axis([0,100,0,0.05]);

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
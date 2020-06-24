
clear
global minf msup

minf=0.01;
msup=10;

m = (0:1e-5:1).*(msup-minf)+minf;

fm05 = fmchab05(m)/integral(@fmrecente, 0.01, 100);
fm03 = fmchab03(m)/integral(@fmchab03, 0.01, 100);
fmkroupa = fm_kroupa(m)/integral(@fm_kroupa, 0.01, 100);
fm_modi = fmchab_modi(m)/integral(@fmchab_modi, 0.01, 100);

figure(1);
set(gca, 'YScale', 'log')
hold on;
plot(log10(m),fmkroupa);
plot(log10(m),fm03);
plot(log10(m),fm05);
plot(log10(m),fm_modi);
hold off

legend('MF1 (Kroupa)','MF2 (2003)','MF3 (2005)', 'fm chabier modif');
ylabel('\xi(log(m))');
xlabel('log_{10} (M_{sol})');
% %axis([0,100,0,0.05]);

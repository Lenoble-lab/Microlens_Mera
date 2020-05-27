
clear

global minf msup

minf=0.01;
msup=5;

m = (0:1e-5:1).*(msup-minf)+minf;

fm05 = fmrecente(m);
fm03 = fmchab03(m);
%fmpower = fm1(m);

figure(1);
hold on;
plot(log10(m),log10(fm05));
plot(log10(m),log10(fm03));
% plot(log10(m),log10(fmpower));

legend('MF1 (2005)','MF2 (2003)');
ylabel('log_{10} (MF)');
xlabel('log_{10} (M)');
% %axis([0,100,0,0.05]);
title('fonction de masse');

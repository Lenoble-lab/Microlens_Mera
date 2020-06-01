clear;
temax = 100;
nbre_bin = temax/5;

tau = 1;gam = 1;te = [1];uT = 1;


exp_ogle_2006
teff_ogle = teff;
exp_macho_2005
teff_macho = teff;
exp_eros_2006
teff_eros = teff;

% length(teff_ogle)
% le = max(
figure(1)
histogram(teff)

%%
rng 'default'

data1 = randn(20,1);
data2 = randn(30,1);
data3 = randn(40,1);
data4 = randn(50,1);

edges = -4:1:4;

h1 = histcounts(data1,edges);
h2 = histcounts(data2,edges);
h3 = histcounts(data3,edges);
h4 = histcounts(data4,edges);

figure
bar(edges(1:end-1),[h1; h2; h3; h4]')

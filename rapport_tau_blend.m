%%Calcul du rapport entre les profondeurs optiques avec et sans blending

%----------------------------------------
% Param�tres de la fonction de luminosité pour le blending (Alibert et al.)
%----------------------------------------
clear; close all;
global Vinf Vsup normfl

Vinf = 22; % magnitudes limites des étoiles observées.
Vsup = 16;

Linf=lum(Vinf);
Lsup=lum(Vsup);

% Norme de la fonction de luminosité :

normfl = integral(@fl, Linf, Lsup);

% Preparation de la table pour le tirage aleatoire des luminosités

ll = (0:1e-5:1).*(Lsup-Linf)+Linf;
fll = probal(ll);
ifll = cumsum(fll);	% primitive 
ifll = ifll-ifll(1);
ifll = real(ifll./ifll(end));	% on fait en sorte que la primitive varie de 0 a 1

f_n_bar = @(n) n.*exp(-n)./(1-exp(-n));

n = 1e5;
n_li = (0:1e-3:1).*20;

%f_li : fraction d'étoiles non blendée
f_li = f_n_bar(n_li);

rapport = zeros(size(n_li));

for k = 1:length(f_li)
    
nbar = n_li (k);
f = f_li(k);

ra = rand(1,n); 
in = find(ra-f> 0);  %Evenements concerné par le blending
out = find(ra-f<= 0);  %évènements non concernés


% Tirage des luminosités (donc des flux)

ra1 = rand(size(in));
ra2 = rand(size(in));

flux1 = interp1(ifll,ll,ra1);
flux2 = interp1(ifll,ll,ra2);

B = flux1 ./ (flux1 + flux2);

% Blending

At=1.34; %on utilse le seuil défini par l'expérience
Umin = rand(size(in));

%Calcul de B_min
Bmin = (At-1)./(ampli(Umin)-1);

i0 = find(B<Bmin);
%Calcul de Uobs
Uobs=maxampli(B.*ampli(Umin)+1-B);

%Calcul du facteur d'amplification
fact(in) = sqrt((maxampli((At - 1)./B +1)).^2-Umin.^2)./(sqrt(1-Uobs.^2));
fact(out) = ones(size(out));
fact(i0)=zeros(size(i0));


il=find(fact~=real(fact)); 
fact(il)=1; % donne parfois des nombres complexes si trop proche de l'amplification infinie, on pose donc une amplification = 1

% On simule te et on observe te,obs, donc nous voulons tracer te,obs

gmean = mean(fact);

gmean_li(k) = gmean;
expo_obs_sur_expo_vrai(k) = nbar/(1-exp(-nbar));
rapport(k) = (gmean * (nbar/(1-exp(-nbar))));

end

figure(1)
hold on;
plot(f_li, rapport)
ylabel('\tau_{obs} / \tau')
xlabel('fraction de sources sans biais de confusion')

figure(2)
plot(f_li, expo_obs_sur_expo_vrai)

figure(3)
plot(f_li, gmean_li)
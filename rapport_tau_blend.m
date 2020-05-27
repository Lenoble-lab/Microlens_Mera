%%Calcul du rapport entre les profondeurs optiques avec et sans blending

%----------------------------------------
% Param�tres de la fonction de luminosité pour le blending (Alibert et al.)
%----------------------------------------

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

n = 1e5;

for k = 1:length(f_li)

ra = rand(n); 
f = f_li(k);
out = find(ra-f> 0);
in = find(ra-f<= 0);

% Tirage des luminosités (donc des flux)

ra1 = rand(size(out));
ra2 = rand(size(in));

flux1 = interp1(ifll,ll,ra1);
flux2 = interp1(ifll,ll,ra2);

B = flux1 ./ (flux1 + flux2);

% Blending

%opération pour pouvoir récupérer les infos plus tard : on conserve les indices et le résultats du blending par rapport à te

blend = zeros(size(in));

At=ampli(uT); %on utilse le seuil défini par l'expérience
Umin = rand(size(in));
Uobs= zeros(size(in));

%Calcul de Uobs
Uobs=maxampli(B.*ampli(Umin)+1-B);

%Calcul du facteur d'amplification
fact = sqrt((maxampli((At - 1)./B +1)).^2-Umin.^2)./(sqrt(1-Uobs.^2));
il=find(fact~=real(fact)); 
fact(il)=1; % donne parfois des nombres complexes si trop proche de l'amplification infinie, on pose donc une amplification = 1
	    % On simule te et on observe te,obs, donc nous voulons tracer te,obs

gmean = mean(fact);

taurblend=taur * gmean * (nbar/(1-exp(-nbar)));
taurblend=real(taurblend);
disp(['tau avec blending (Alibert 2005)  = ' num2str(taurblend)]);

end
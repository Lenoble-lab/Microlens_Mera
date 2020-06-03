%----------------------------------------------------
%Nécessite en argument nbar et f (fraction d'évènements non blendé)
%----------------------------------------------------

% Preparation de la table pour le tirage aleatoire des luminosités


ll = (0:1e-5:1).*(Lsup-Linf)+Linf;
fll = probal(ll);
ifll = cumsum(fll);	% primitive 
ifll = ifll-ifll(1);
ifll = real(ifll./ifll(end));	% on fait en sorte que la primitive varie de 0 a 1

% Tirage est évenements concernés

ra = rand(size(te)); 
in = find(ra-f> 0);  %Evenements concerné par le blending
out = find(ra-f<= 0);  %évènements non concernés


% Tirage des luminosités (donc des flux)

ra1 = rand(size(te(in)));
ra2 = rand(size(te(in)));

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

%Calcul de B_min
Bmin = (At-1)./(ampli(Umin)-1);

i0 = find(B<Bmin);
fact(i0)=zeros(size(i0));

%Récupération des résultats
teblend = te;
teblend(in)=te(in).*fact; % on a appliqué le blending à te et on a conservé l'ordre de te (important pour le blending après efficacité)

%Pour calculer gmean
fact(out) = ones(size(out));

gmean = mean(fact);

taurblend=taur * gmean * (nbar/(1-exp(-nbar)));
taurblend=real(taurblend);
disp(['tau avec blending (Alibert 2005)  = ' num2str(taurblend)]);

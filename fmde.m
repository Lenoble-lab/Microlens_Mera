% Fonction de masse dN/dm des lentilles.

function pm = fmde(m)
global minfde msupde 
pm = zeros(size(m));

global minf msup

%---
%MF1
%---

A=0.053;
mbar=1.1346;
beta=0.6;
alpha=2.7;


%i1 = find(m>=minf & m<=msup);
%if (length(i1)>=1)
%  pm(i1) = A.*exp(-(mbar./m(i1)).^beta).*(m(i1).^(-alpha));
%end



%-------
%MF1 bis
%-------

A=0.126;
mbar=5.49;
beta=0.39;
alpha=2.7;


%i1 = find(m>=minf & m<=msup);
%if (length(i1)>=1)
%  pm(i1) = A.*exp(-(mbar./m(i1)).^beta).*(m(i1).^(-alpha));
%end



%--------
%MF1 tris
%--------

A=0.694;
mbar=76.62;
beta=0.3;
alpha=3.2;


%i1 = find(m>=minf & m<=msup);
%if (length(i1)>=1)
%  pm(i1) = A.*exp(-(mbar./m(i1)).^beta).*(m(i1).^(-alpha));
%end


%----------
%MF1 quater
%----------

A=3.0;
mbar=716.4;
beta=0.25;
alpha=3.3;


%i1 = find(m>=minf & m<=msup);
%if (length(i1)>=1)
%  pm(i1) = A.*exp(-(mbar./m(i1)).^beta).*(m(i1).^(-alpha));
%end




%---
%MF2
%---

alpha=-1.24;

%i1 = find(m>=minf);
%if (length(i1)>=1)
%  pm(i1) = m(i1).^(alpha);
%end


%-------
%MF2 bis
%-------

alpha=-1.55;

%i1 = find(m>=minf);
%if (length(i1)>=1)
%  pm(i1) = m(i1).^(alpha);
%end


%--------
%MF2 tris
%--------

alpha=-1.5;

%i1 = find(m>=minf);
%if (length(i1)>=1)
%  pm(i1) = m(i1).^(alpha);
%end


%---
%MF3
%---

sigma=0.677;
m0=0.095;

%i1 = find(m>=minf);
%if (length(i1)>=1)
%  pm(i1) =exp(-((log10(m)-log10(m0)).^2)./(2*sigma^2))./m;
%end  

%---------------------------
%rajout d'une quantite de WD 
%---------------------------

%largeur=0.01;
%massWD=0.6;
%fraction=0;

%i=find((m<massWD+0.5*largeur) & (m>massWD-0.5*largeur));
%if(length(i)>=1)
%  pm(i)=pm(i)+fraction;
%end  


%---------------------------
%fonctions de masse externes
%---------------------------

pm=fmrecente(m);
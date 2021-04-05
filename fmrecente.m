% Fonction de masse dN/dm des lentilles.
%C'est cette fonction qu'on utilise dans le code main, d'où l'intéret de
%changer uniquement la fonction qu'elle renvoie. On peut aussi aller
%directement dans les fonctions de masse du halo/disque/bulbe si on veut
%différencier les fonction selon les populations considérées
function pm = fmrecente(m)

% pm = fmchab03(m);
pm = fmchab05(m);
% pm = fm_kroupa_modif(m);
% pm = fmchab_modi(m);
% pm = fm_basu_rana(m);
% pm = fm_kroupa(m);
% pm = fmchab(m);
% pm = fmchab_2014_cas_1(m);

% pm = zeros(size(m));
% % [~,idx] = find(abs(m-1)<0.3);
% [~,idx] = min(abs(m-1));
% pm(idx) = 1;
end

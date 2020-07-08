% Fonction de masse dN/dm des lentilles.
% cette fonction de masse est la plus recente : 2005

function pm = fmrecente(m)

% pm = fmchab03(m);
pm = fmchab05(m);
% pm = fm_kroupa_modif(m);
% pm = fmchab_modi(m);
% pm = PDMF_Maraston(pm,m);
% pm = PDMF_gould(pm,m);

% pm = zeros(size(m));
% % [~,idx] = find(abs(m-1)<0.3);
% [~,idx] = min(abs(m-1));
% pm(idx) = 1;
end

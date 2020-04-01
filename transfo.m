function res = transfo ( t )  

% sous-fonction qui permet de transformer en tableau ayant plusieurs meme
% elements en un tableau ou tout les elements n'apparaissent qu'une seule
% fois (c'est pour la fonction Matlab interp1)

[m,l] = size(t);

if (m>1)
    disp(['erreur'])'
else
    if ( l == 0 )
        res = [];
    elseif ( l == 1 )
        res = [t(1)];
    else
        res = [t(1)];
        for i = 2:l;
            if ( ~(appartient(res,t(i))) )
                res=[res,t(i)];
            end;
        end;
    end;
end;

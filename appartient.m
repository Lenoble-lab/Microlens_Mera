function res = appartient (t,a)

% fonction qui teste si l'element a appartient au tableau t et renvoie i si
% oui et 0 si non

l=length(t);
i=1;
while ( i<=l & t(i)~=a )
    i=i+1;
end;
if ( i>l )
    res=0;
else
    res=1;
end;


function [tt,effic] = tri ( t , eff )

% fonction qui trie le tableau t par elements croissant et effectue la
% meme operation sur le tableau eff

l = length(t);

tt=t;
effic=eff;

for i=1:l;
    for j=i:l;
        if ( tt(j) < tt(i) )
            c=tt(i);
            tt(i)=tt(j);
            tt(j)=c;
            c=effic(i);
            effic(i)=effic(j);
            effic(j)=c;
        end;
    end;
end;





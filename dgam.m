function res = dgam(x,L,v,m)

global cosb cosbl sinb sinl siglr siglt siglp L0 Ro vlimit
% mise en forme des inputs
i = find(x<0 | x>1 | v < - vlimit | v > vlimit);

px = rholens(x.*L).*sqrt(x.*(1-x));	% terme rho(xL)*sqrt(x(1-x))

res = px.*v.*sqrt(m).*L.^1.5;	% Calcul du dgamma (cf eq 3.25 p. 57 de ma these)

%pas de terme nsource(L) car les distances
%sont tirees au hasard avec une loi proportionnelle
%a nsource



if (length(i)>=1) res(i) = zeros(1,length(i));end


    
function res = dgamma(x,glr,glt,glp,L,gsr,gst,gsp)

global cosb cosbl sinb sinl siglr siglt siglp L0 Ro
% mise en forme des inputs
i = find(x<0 | x>1 | glr>=0.99999 | glr<=0.00001 | glp>=0.99999 | glp<=0.00001 | glt>=0.99999 | glt<=0.00001 | gsr>=0.99999 | gsr<=0.00001 | gsp>=0.99999 | gsp<=0.00001 | gst>=0.99999 | gst<=0.00001  );

[R, z, th] = toGC(x.*L); 
px = rholens(x.*L).*sqrt(x.*(1-x));	% terme rho(xL)*sqrt(x(1-x))

v = vperp(x,glr,glt,glp,L,gsr,gst,gsp);
res = px.*v.*L.^1.5;	% Calcul du dgamma (cf eq 3.25 p. 57 these Alibert)
if (length(i)>=1) res(i) = zeros(1,length(i));end


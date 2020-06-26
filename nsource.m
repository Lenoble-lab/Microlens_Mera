function res = nsource (x)

global dsup dinf 

res = zeros(size(x));

bet=0;


i1 = find(x>=dinf & x<=dsup);
if (length(i1)>=1)
  res (i1) = denssource(x(i1)).*x(i1).^(2.*bet+2);
end
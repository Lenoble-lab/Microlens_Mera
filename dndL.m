function res = dndL (x)

global beta Lmax Lmin

res = zeros(size(x));

i1 = find(x>=Lmin & x<=Lmax);
if (length(i1)>=1)
  res (i1) = x(i1).^(2*beta);
end



function semilogxhist(val,M)
% semilogxhist - generate histogram with M bars and log-scale x axis
vmin=min(val); vmax=max(val);
edges=vmin*(vmax/vmin).^([0:M]/M);
count=histc(val,edges); 
if size(count,2)==1, count=count'; end 
x=edges(sort([1:M 1:M])); 
y=[0 count(sort([1:M-1 1:M-1])) 0];

% outline only: semilogx(x, y, '-');
plot(x, y, '-'); 
% fill(x, y, 'b'); 
set(gca,'XScale','log');
end
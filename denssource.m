function res = denssource(x)

global mmeande mmeandm mmeanbu mmeanh 

[R, z, th] = toGC(x);

res=zeros(size(R));

%---------------
% source : bulbe
%---------------

res = res + rhobulbe(R,z,th)*mmeanbu;

%----------------------
% source : disque mince
%----------------------

res = res + rhodm(R,z,th)*mmeandm;

%----------------------
% source : disque epais
%----------------------

res = res + rhode(R,z,th)*mmeande;

%--------------
% source : halo
%--------------

%res = res + rhohalo(R,z,th)*mmeanh;

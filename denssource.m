function res = denssource(x)

[R, z, th] = toGC(x);

res=zeros(size(R));

%---------------
% source : bulbe
%---------------

% res = res + rhobulbe(R,z,th)*mmeanbu;
res = res + rhobulbe(R,z,th);

%----------------------
% source : disque HetG
%----------------------

% res = res + rhodHetG(R, z, th);

%----------------------
% source : disque mince
%----------------------

% res = res + rhodm(R,z,th)*mmeandm;
res = res + rhodm(R,z,th);

%----------------------
% source : disque epais
%----------------------

% res = res + rhode(R,z,th)*mmeande;
% res = res + rhode(R,z,th);

%--------------
% source : halo
%--------------

%res = res + rhohalo(R,z,th)*mmeanh;

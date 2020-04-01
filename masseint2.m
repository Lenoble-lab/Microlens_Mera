function int  = masseint2 ()
% calcul d'integral triple

rmax = 1e5;
zmax = 1e5;
x=2*pi;
pr=1e0;
pz=1e0;
pt=1e0;

%int = triplequad('rhozhaoint',0,1e20,-1e20,1e20,-x,x);

%int = quadl( ( quadl( ( quadl('rhozhaoint',0,x,pt) ) ,-zmax,zmax,pz) ) ,0,rmax,pr);

R = 0:pr:rmax;
th = 0:pt:x;
z = -zmax:pz:zmax;

lr = length(R);
lt = length(th);
lz = length(z);

int=0;
for i = 1:lr;
    for j = 1:lt;
        for k = 1:lz;
            int = int + rhozhaoint(R(i),z(k),th(j));
        end;
    end;
end;

int = int*pr*pt*pz

disp(num2str(int));












% 5.0326e+08
% 5.4669e+07

% 8.4227e+07
% 9.2429e+07
% 1.0538e+08
% 1.2443e+08
% 1.2608e+08
% 15396331.8563



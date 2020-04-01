function int  = masse ()
x=pi;
%int = triplequad('rhozhao',0.000000001,100000000,-100000000,100000000,0,x)
int = triplequad('rhozhao',1e-20,1e20,-1e20,1e20,-x,x);
disp(num2str(int));

% 5.0326e+08
% 5.4669e+07

% 8.4227e+07
% 9.2429e+07
% 1.0538e+08
% 1.2443e+08
% 1.2608e+08
% 15396331.8563


%42132910.4981
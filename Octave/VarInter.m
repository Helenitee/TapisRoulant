function [F_c,omega,tcF,tcK]=VarInter(k,l_0,m,g,v,nu)
%% Variables intermediares
% qui dependent des variables d'entrees

F_c = nu*m*g; % max Force entrainement [N]

omega = sqrt(k/m); % pulsation         [s^-1]

tcF = F_c/(k*v);

tcK = 2*pi/omega;

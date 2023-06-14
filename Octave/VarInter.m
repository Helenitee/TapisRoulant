function [F_c,w,tcF,tcK]=VarInter(k,l_0,m,g,v,nu)
  F_c = nu*m*g;      % max Force entrainement                  [N]
  w = sqrt(k/m); % pulsation                               [s^-1]
  tcF = F_c/(k*v);   % Temps carct deter par condition ini     [s]
  tcK = 2*pi/w;  % Temps carct deter par pusation          [s]
%endfunction

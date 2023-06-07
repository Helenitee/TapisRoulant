% dérivée de x en régime de Glissement
function[dot_x] = d_xG(t,t_d,x_d,v_d,omega,phi)
  dot_x = - (x_d - phi) * omega * sin(omega * (t - t_d)) + v_d * cos(omega * (t - t_d));
endfunction

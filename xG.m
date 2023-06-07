% x en régime de Glissement
function x = xG(t,t_d,x_d,v_d,v,omega,phi)
  x = (x_d - phi) * cos(omega * (t - t_d)) + v_d/omega * sin(omega * (t - t_d)) + phi;
endfunction

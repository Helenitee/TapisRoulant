% Vitesse ; reg Glissement
function dot_x = vG(t,t_d,x_d,v_d,w,phi)
  dot_x = - (x_d - phi) * w * sin(w * (t - t_d)) + v_d * cos(w * (t - t_d));
endfunction

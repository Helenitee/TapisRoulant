% Position ; reg Glissement
function x = xG(t,t_d,x_d,v_d,w,phi)
x = (x_d - phi) * cos(w * (t - t_d)) + v_d/w * sin(w * (t - t_d)) + phi;
endfunction

% Position en regime d'adherence
function x = xA(t,t_d,x_d,v)
  x = v * (t - t_d) + x_d;
end

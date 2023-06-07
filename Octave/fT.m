function f = fT(t,t_d,x_d,v_d,v,k,l_0,F_c,re)
  if(v_d ~= v)
    s = - sign(v_d - v);
  else
    s = 1;
  end
  %
  if(re == 'ad')
    x = xA(t,t_d,x_d,v);
    f = k * (x - l_0);
  elseif(re == 'gl')
    f = s * F_c;
  end
end

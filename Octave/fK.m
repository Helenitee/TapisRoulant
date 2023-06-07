function f = fK(t,t_d,x_d,v_d,v,k,l_0,omega,phi,re)
  if(re=='ad')
    x = xA(t,t_d,x_d,v);
  elseif(re=='gl')
    x=xG(t,t_d,x_d,v,v_d,omega,phi);
  end
   f = - k * (x - l_0);
endfunction

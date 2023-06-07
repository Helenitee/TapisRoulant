function phi = Phi(v_d,k,l_0,v,F_c)
  if(v_d ~= v)
    s = sign(v_d - v);
  else
    s = 1;
  end
  phi = l_0 + s * F_c/k;
endfunction

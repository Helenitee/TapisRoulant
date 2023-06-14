function lst = List(min_x,max_x,min_v,max_v,nx,nv,k,l_0,m,g,v,nu,F_c,w,tcF,tcK)
  lst = [];
  nt = 200;
  t_0 = 0;
  
  x0 = linspace(min_x,max_x,nx);
  v0 = linspace(min_v,max_v,nv);

  for i = 1:nx;
    x_0 = x0(i);
    
    for j = 1:nv;
      v_0 = v0(j);

      %%Phase 1 : glissement ----------------------
      phi = Phi(v_0,k,l_0,v,F_c);
      t01=linspace(t_0,t_0 + 4*max(tcF,tcK),nt);
      
      CostG = @(t,t_d,x_d,v_d,w,phi,v) ((vG(t,t_d,x_d,v_d,w,phi) - v).^2);
      it1=find(diff(sign(diff(CostG(t01,t_0,x_0,v_0,w,phi,v))))==2,1);

      t_1 = fminsearch(@(t) CostG(t,t_0,x_0,v_0,w,phi,v),t01(it1+1));
      x_1 = xG(t_1,t_0,x_0,v_0,w,phi);
      
      %% Phase 2 : adherence ------------------------
      v_1 = v;
      phi = Phi(v_1,k,l_0,v,F_c);

      t12 = linspace(t_1,t_1 + 4*max(tcK,tcF),nt);
      
      CostA = @(t,t_d,x_d,v_d,v,k,l_0,F_c) (abs(fT(t,t_d,x_d,v_d,v,k,l_0,F_c,'ad') - F_c).^2);
      it2 = find(diff(sign(diff(CostA(t12,t_1,x_1,v_1,v,k,l_0,F_c))))==2,1);
      
      t_2 = fminsearch(@(t) CostA(t,t_1,x_1,v_1,v,k,l_0,F_c),t12(it2+1));
      x_2 = xA(t_2,t_1,x_1,v_1);
      
      %% Phase 3 : glissement ------------------------
      v_2 = v;
      phi = Phi(v_2,k,l_0,v,F_c);

      t23 = linspace(t_2,t_2 + 2*max(tcF,tcK),nt);
      
      it3 = find(diff(sign(diff(CostG(t23,t_2,x_2,v_2,w,phi,v))))==2,1);
      
      t_3 = fminsearch(@(t) CostG(t,t_2,x_2,v_2,w,phi,v),t23(it3+1));
      x_3 = xG(t_3,t_2,x_2,v_2,w,phi);

      lst = [lst;x_0,v_0,t_1,t_2,t_3,t_2-t_1];%X_1,x_2,x_3];
    end
  end
endfunction

  close all; clear all; clc;

 %% Variables d'entree
[k,l_0,m,g,v,nu]=VarEntree;

%% Variables intermediaires
[F_c,w,tcF,tcK]=VarInter(k,l_0,m,g,v,nu);

%% CAS 2.1 :
 lst = List(l_0+0.11,l_0+0.24,v,v,14,1,k,l_0,m,g,v,nu,F_c,w,tcF,tcK);
 x1 = lst(:,1);
 y1 = lst(:,3);
 invy1 = 1./y1;
 p1 = polyfit(x1,invy1,4);
 reg1 = polyval(p1,x1);

 figure(1)
 hold on
 plot(x1,1./reg1,'-','color',[0 0 .5],'LineWidth',1)
 plot(x1,y1,'+','color',[1 0 .5],'LineWidth',.7)
 grid minor;
 legend('modèle','valeurs calculées','Interpreter', 'latex','location', 'northeast','fontsize', 20);
 xlabel('v_0');
 ylabel('t_1');
 title('Représentation de t_1 suivant v_0 ', 'fontsize',25);

 x2 = lst(:,1);
 y2 = lst(:,6);

 p2 = polyfit(x2,y2,1);
 reg2 = polyval(p2,x2);


 figure(2)
 hold on
 plot(x2,reg2,'-','color',[0 0 .5],'LineWidth',1)
 plot(x2,y2,'+','color',[1 0 .5],'LineWidth',.7)
 grid minor;
 legend('modèle','valeurs calculées','location', 'northeast','fontsize', 20);
 xlabel('v_0');
 ylabel('t_2 - t_1');
 title('Représentation de t_2 - t_1 suivant v_0 ', 'fontsize', 25);

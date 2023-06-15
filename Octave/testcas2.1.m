 close all; clear all; clc;

 %% Variables d'entree
[k,l_0,m,g,v,nu]=VarEntree;

%% Variables intermediaires
[F_c,w,tcF,tcK]=VarInter(k,l_0,m,g,v,nu);

%% CAS 2.1 :
 lst = List(l_0,l_0,-1,1,1,21,k,l_0,m,g,v,nu,F_c,w,tcF,tcK);
 x11 = lst(1:11,2);
 x12 = lst(13:21,2);
 y11 = lst(1:11,3);
 y12 = lst(13:21,3);
 x1 = [x11;x12];
 y1 = [y11;y12];
 p11 = polyfit(x11,y11,2);
 reg11 = polyval(p11,x11);
 p12 = polyfit(x12,y12,2);
 reg12 = polyval(p12,x12);

 figure(1)
 hold on
 plot(x11,reg11,'-','color',[0 0 .5],'LineWidth',1)
 plot(x12,reg12,'-','color',[0 .5 .7],'LineWidth',1)
 plot(x1,y1,'+','color',[1 0 .5],'LineWidth',.7)
 grid minor;
 legend('-0.0415 x^2 - 0.1047 x + 0.0103','- 0.0405 x^2 + 0.1136 x - 0.0110','valeurs calculées','Interpreter', 'latex','location', 'northeast','fontsize', 20);
 xlabel('v_0');
 ylabel('t_1');
 title('Représentation de t_1 suivant v_0 ', 'fontsize',25);

 x21 = lst(1:11,2);
 x22 = lst(13:21,2);
 y21 = lst(1:11,6);
 y22 = lst(13:21,6);
 x2 = [x21;x22];
 y2 = [y21;y22];
 p21 = polyfit(x21,y21,2);
 reg21 = polyval(p21,x21);
 p22 = polyfit(x22,y22,2);
 reg22 = polyval(p22,x22);

 figure(2)
 hold on
 plot(x21,reg21,'-','color',[0 0 .5],'LineWidth',1)
 plot(x22,reg22,'-','color',[0 .5 .7],'LineWidth',1)
 plot(x2,y2,'+','color',[1 0 .5],'LineWidth',.7)
 grid minor;
 legend('0.2844 x^2 - 0.1003 x + 0.4789','- 0.2315 x^2 - 0.1690 x + 0.5213','valeurs calculées','location', 'northeast','fontsize', 20);
 xlabel('v_0');
 ylabel('t_2 - t_1');
 title('Représentation de t_2 - t_1 suivant v_0 ', 'fontsize', 25);

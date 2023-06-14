
close all; clear all; clc;
format short

%% Variables d'entree
[k,l_0,m,g,v,nu]=VarEntree;

%% Variables intermediaires
[F_c,w,tcF,tcK]=VarInter(k,l_0,m,g,v,nu);

%% Variables initialisation
t_0 = 0   % tmps ini  [s]

% Cas 2 -------------------------------------
x_0 = l_0 + 0.1;
v_0 = v;
% Cas 3 -------------------------------------
%x_0 = l_0 + 0.1;
%v_0 = v;
% Cas 4 -------------------------------------
%x_0 = l_0 + 0.1;
%v_0 = v + 0.2;
% ------------------------------------------
nt=200;

%%Phase 1 : glissement ----------------------
phi = Phi(v_0,k,l_0,v,F_c);

t1=t_0+4*max(tcF,tcK);
t01 = linspace(t_0,t1,nt);
x01 = xG(t01,t_0,x_0,v_0,w,phi);
v01 = vG(t01,t_0,x_0,v_0,w,phi);
fk01 = fK(t01,t_0,x_0,v_0,v,k,l_0,w,phi,'gl');
ft01 = fT(t01,t_0,x_0,v_0,v,k,l_0,F_c,'gl');

figure(1);
subplot(3,1,1);hold on;
plot(t01,l_0*ones(size(t01)),'-k')
plot(t01,(l_0+F_c/k)*ones(size(t01)),'--k')
plot(t01,(l_0-F_c/k)*ones(size(t01)),'--k')
plot(t01,x01,'-b')
subplot(3,1,2);hold on;
plot(t01,v*ones(size(t01)),'-k')
plot(t01,v01,'-b')
subplot(3,1,3);hold on;
plot(t01,-F_c*ones(size(t01)),'-k')
plot(t01,+F_c*ones(size(t01)),'-k')
plot(t01,fk01,'-b')
plot(t01,ft01,'-r')

%% 
CostG = @(t,t_0,x_0,v_0,w,phi,v) ((vG(t,t_0,x_0,v_0,w,phi) - v).^2);
figure(11);
subplot(2,1,1);hold on;
plot(t01,CostG(t01,t_0,x_0,v_0,w,phi,v),'b')
subplot(2,1,2);hold on;
stem(t01(1:nt-2),diff(sign(diff(CostG(t01,t_0,x_0,v_0,w,phi,v)))),'xb')
it1=find(diff(sign(diff(CostG(t01,t_0,x_0,v_0,w,phi,v))))==2,1);
subplot(2,1,1);hold on;
plot(t01(it1+1),CostG(t01(it1+1),t_0,x_0,v_0,w,phi,v),'o')

%% 
t_1 = fminsearch(@(t) CostG(t,t_0,x_0,v_0,w,phi,v),t01(it1+1))
t_01 = linspace(t_0,t_1,nt);
x_01 = xG(t_01,t_0,x_0,v_0,w,phi);
v_01 = vG(t_01,t_0,x_0,v_0,w,phi);
fk_01 = fK(t_01,t_0,x_0,v_0,v,k,l_0,w,phi,'gl');
ft_01 = fT(t_01,t_0,x_0,v_0,v,k,l_0,F_c,'gl');
figure(1);
subplot(3,1,1);hold on;
plot(t_01,l_0*ones(size(t_01)),'-k')
plot(t_01,(l_0+F_c/k)*ones(size(t_01)),'--k')
plot(t_01,(l_0-F_c/k)*ones(size(t_01)),'--k')
plot(t_01,x_01,'-b','LineWidth',4)
subplot(3,1,2);hold on;
plot(t_01,v*ones(size(t_01)),'-k')
plot(t_01,v_01,'-b','LineWidth',4)
subplot(3,1,3);hold on;
plot(t_01,-F_c*ones(size(t_01)),'-k')
plot(t_01,+F_c*ones(size(t_01)),'-k')
plot(t_01,fk_01,'-b','LineWidth',4)
plot(t_01,ft_01,'-r','LineWidth',4)

%%
figure(100);
subplot(3,1,1);hold on;
plot(t_01,l_0*ones(size(t_01)),'-k')
plot(t_01,(l_0+F_c/k)*ones(size(t_01)),'--k')
plot(t_01,(l_0-F_c/k)*ones(size(t_01)),'--k')
plot(t_01,x_01,'-b','LineWidth',4)
subplot(3,1,2);hold on;
plot(t_01,v*ones(size(t_01)),'-k')
plot(t_01,v_01,'-b','LineWidth',4)
subplot(3,1,3);hold on;
plot(t_01,-F_c*ones(size(t_01)),'-k')
plot(t_01,+F_c*ones(size(t_01)),'-k')
plot(t_01,fk_01,'-b','LineWidth',4)
plot(t_01,ft_01,'-r','LineWidth',4)

%%Phase 2 : adherence----------------------
x_1=xG(t_1,t_0,x_0,v_0,w,phi);
v_1=v;
phi = Phi(v_1,k,l_0,v,F_c);

t2=t_1+4*max(tcF,tcK);
t12 = linspace(t_1,t2,nt);
x12 = xA(t12,t_1,x_1,v_1);
v12 = v*ones(size(t12));
fk12 = fK(t12,t_1,x_1,v_1,v,k,l_0,w,phi,'ad');
ft12 = fT(t12,t_1,x_1,v_1,v,k,l_0,F_c,'ad');

%%
figure(2);
subplot(3,1,1);hold on;
plot(t12,l_0*ones(size(t12)),'-k')
plot(t12,(l_0+F_c/k)*ones(size(t12)),'--k')
plot(t12,(l_0-F_c/k)*ones(size(t12)),'--k')
plot(t12,x12,'-b')
subplot(3,1,2);hold on;
plot(t12,v*ones(size(t12)),'-k')
plot(t12,v*ones(size(t12)),'-k')
plot(t12,v12,'-b')
subplot(3,1,3);hold on;
plot(t12,-F_c*ones(size(t12)),'-k')
plot(t12,+F_c*ones(size(t12)),'-k')
plot(t12,-F_c*ones(size(t12)),'-k')
plot(t12,+F_c*ones(size(t12)),'-k')
plot(t12,fk12,'-b')
plot(t12,ft12,'-r')

%% 
CostG = @(t,t_1,x_1,v_1,v,k,l_0,F_c) (abs(fT(t,t_1,x_1,v_1,v,k,l_0,F_c,'ad') - F_c).^2);
figure(12);
subplot(2,1,1);hold on;
plot(t12,CostG(t12,t_1,x_1,v_1,v,k,l_0,F_c),'b')
subplot(2,1,2);hold on;
stem(t12(1:nt-2),diff(sign(diff(CostG(t12,t_1,x_1,v_1,v,k,l_0,F_c)))),'xb')
it2=find(diff(sign(diff(CostG(t12,t_1,x_1,v_1,v,k,l_0,F_c))))==2,1);
subplot(2,1,1);hold on;
plot(t12(it2+1),CostG(t12(it2+1),t_1,x_1,v_1,v,k,l_0,F_c),'o')

%% 
t_2 = fminsearch(@(t) CostG(t,t_1,x_1,v_1,v,k,l_0,F_c),t12(it2+1))
t_12 = linspace(t_1,t_2,nt);
x_12 = xA(t_12,t_1,x_1,v_1);
v_12 = v*ones(size(t_12));
fk_12 = fK(t_12,t_1,x_1,v_1,v,k,l_0,w,phi,'ad');
ft_12 = fT(t_12,t_1,x_1,v_1,v,k,l_0,F_c,'ad');
figure(2);
subplot(3,1,1);hold on;
plot(t_12,x_12,'-b','LineWidth',4)
subplot(3,1,2);hold on;
plot(t_12,v_12,'-b','LineWidth',4)
subplot(3,1,3);hold on;
plot(t_12,fk_12,'-b','LineWidth',4)
plot(t_12,ft_12,'-r','LineWidth',4)

%
figure(100);
subplot(3,1,1);hold on;
plot(t_12,l_0*ones(size(t_12)),'-k')
plot(t_12,(l_0+F_c/k)*ones(size(t_12)),'--k')
plot(t_12,(l_0-F_c/k)*ones(size(t_12)),'--k')
plot(t_12,x_12,'-b','LineWidth',4)
subplot(3,1,2);hold on;
plot(t_12,v*ones(size(t_12)),'-k')
plot(t_12,v_12,'-b','LineWidth',4)
subplot(3,1,3);hold on;
plot(t_12,-F_c*ones(size(t_12)),'-k')
plot(t_12,+F_c*ones(size(t_12)),'-k')
plot(t_12,fk_12,'-b','LineWidth',4)
plot(t_12,ft_12,'-r','LineWidth',4)

%%Phase 1 : glissement ----------------------
x_2 = xA(t_2,t_1,x_1,v_1);
v_2 = v;
phi = Phi(v_2,k,l_0,v,F_c);

t3=t_2+2*max(tcF,tcK);
t23 = linspace(t_2,t3,nt);
x23 = xG(t23,t_2,x_2,v_2,w,phi);
v23 = vG(t23,t_2,x_2,v_2,w,phi);
fk23 = fK(t23,t_2,x_2,v_2,v,k,l_0,w,phi,'gl');
ft23 = fT(t23,t_2,x_2,v_2,v,k,l_0,F_c,'gl');

figure(3);
subplot(3,1,1);hold on;
plot(t23,l_0*ones(size(t23)),'-k')
plot(t23,(l_0+F_c/k)*ones(size(t23)),'--k')
plot(t23,(l_0-F_c/k)*ones(size(t23)),'--k')
plot(t23,x23,'-b')
subplot(3,1,2);hold on;
plot(t23,v*ones(size(t23)),'-k')
plot(t23,v23,'-b')
subplot(3,1,3);hold on;
plot(t23,-F_c*ones(size(t23)),'-k')
plot(t23,+F_c*ones(size(t23)),'-k')
plot(t23,fk23,'-b')
plot(t23,ft23,'-r')

figure(100);
subplot(3,1,1);hold on;
plot(t23,l_0*ones(size(t23)),'-k')
plot(t23,(l_0+F_c/k)*ones(size(t23)),'--k')
plot(t23,(l_0-F_c/k)*ones(size(t23)),'--k')
plot(t23,x23,'-b','LineWidth',4)
subplot(3,1,2);hold on;
plot(t23,v*ones(size(t23)),'-k')
plot(t23,v23,'-b','LineWidth',4)
subplot(3,1,3);hold on;
plot(t23,-F_c*ones(size(t23)),'-k')
plot(t23,+F_c*ones(size(t23)),'-k')
plot(t23,fk23,'-b','LineWidth',4)
plot(t23,ft23,'-r','LineWidth',4)
% % -------------------------------------------
% % Temps
% %t_1 = fminsearch(@(t) Cost1(t), .2*tcK)       % Cas 2
% %t_1 = fminsearch(@(t) Cost1(t),t_0+1/2*tcK)  % Cas 3
% %t_1 = fminsearch(@(t) Cost1(t),t_0+tcK)      % Cas 4
% t_01 = linspace(t_0,t_1,nt);
% t01=linspace(t_0,t_0+4*max(tcF,tcK),nt);
% % Position
% x_1 = xG(t_1,t_0,x_0,v_0,v,w,phi_1)
% x_01 = xG(t_01,t_0,x_0,v_0,v,w,phi_1);
% x01 = xG(t01,t_0,x_0,v_0,v,w,phi_1);
% % Vitesse
% v_1 = d_xG(t_1,t_0,x_0,v_0,w,phi_1);
% v_01 = d_xG(t_01,t_0,x_0,v_0,w,phi_1);
% v01 = d_xG(t01,t_0,x_0,v_0,w,phi_1);
% % Force ressort
% f_1k = fK(t_1,t_0,x_0,v_0,v,k,l_0,w,phi_1,'gl');
% f_01k = fK(t_01,t_0,x_0,v_0,v,k,l_0,w,phi_1,'gl');
% % Force tapis
% f_1t = fT(t_1,t_0,x_0,v_0,v,k,l_0,F_c,'gl');
% f_01t = fT(t_0,t_0,x_0,v_0,v,k,l_0,F_c,'gl') * ones(size(t_01));
% 
% 
% %% Phase 2 : adherence ------------------------
% phi = Phi(v,k,l_0,v,F_c);
% Cost2=@(t) (abs(fT(t,t_1,x_1,v,v,k,l_0,F_c,'ad')) - F_c).^2;
% % ---------------------------------------------
% t_2 = t_1+tcF-(x_1-l_0)/v
% t_12 = linspace(t_1,t_2,nt);
% t12=linspace(t_1,t_1 + 4*max(tcF,tcK),nt);
% % Position
% x_2 = xA(t_2,t_1,x_1,v)
% x_12 = xA(t_12,t_1,x_1,v);
% x12 = xA(t12,t_1,x_1,v);
% % Vitesse
% v_12 = v*ones(size(t_12));
% v12 = v*ones(size(t12));
% % Force du tapis
% f_2t =  fT(t_2,t_1,x_1,v,v,k,l_0,F_c,'ad');
% f_12t = fT(t_12,t_1,x_1,v,v,k,l_0,F_c,'ad');
% % Force de rappel :
% f_2k = fK(t_2,t_1,x_1,v,v,k,l_0,w,phi,'ad');
% f_12k = fK(t_12,t_1,x_1,v,v,k,l_0,w,phi,'ad');
% 
% 
% %% Phase 3 : glissement ------------------------
% Cost3 = @(t) ((d_xG(t,t_2,x_2,v,w,phi) - v).^2);
% % -------------------------------------------
% t_3 = fminsearch(@(t) Cost3(t),t_2 + tcK)    % Cas2
% %t_3 = fminsearch(@(t) Cost3(t),t_2 + tcK)   % Cas3
% %t_3 = fminsearch(@(t) Cost3(t),t_2 + tcK)   % Cas4
% t_23 = linspace(t_2,t_3,nt);
% t23=linspace(t_2,t_2+4*max(tcF,tcK),nt);
% % Position
% x_3 = xG(t_3,t_2,x_2,v,v,w,phi);
% x_23 = xG(t_23,t_2,x_2,v,v,w,phi);
% x23 = xG(t23,t_2,x_2,v,v,w,phi);
% % Vitesse
% v_3 = d_xG(t_3,t_2,x_2,v,w,phi);
% v_23 = d_xG(t_23,t_2,x_2,v,w,phi);
% v23 = d_xG(t23,t_2,x_2,v,w,phi);
% % Force ressort
% f23k = fK(t23,t_2,x_2,v,v,k,l_0,w,phi,'gl');
% f_23k = fK(t_23,t_2,x_2,v,v,k,l_0,w,phi,'gl');
% f_3k = fK(t_3,t_2,x_2,v,v,k,l_0,w,phi,'gl');
% % Force tapis
% f23t = fT(t23,t_2,x_2,v,v,k,l_0,F_c,'gl') * ones(size(t23));
% f_23t = fT(t_23,t_2,x_2,v,v,k,l_0,F_c,'gl') * ones(size(t_23));
% f_3t = fT(t_3,t_2,x_2,v,v,k,l_0,F_c,'gl');
% 
% 
% %% Phase 4 : adherence ---------------------
% Cost4=@(t) (abs(fT(t,t_3,x_3,v,v,k,l_0,F_c,'ad')) - F_c).^2;
% % -------------------------------------------
% t_4 = t_3+tcF - (x_3 - x_0)/v
% t_34 = linspace(t_3,t_4,nt);
% t34=linspace(t_3,t_3 + 4*max(tcF,tcK),nt);
% % Position
% x_4 = xA(t_4,t_3,x_3,v);
% x_34 = xA(t_34,t_3,x_3,v);
% x34 = xA(t34,t_3,x_3,v);
% % Vitesse
% v_34 = v*ones(size(t_34));
% v34 = v*ones(size(t34));
% % Force du tapis
% f_34t = fT(t_34,t_3,x_3,v,v,k,l_0,F_c,'ad');
% f_4t = fT(t_4,t_3,x_3,v,v,k,l_0,F_c,'ad');
% % Force de rappel :
% f_34k = fK(t_34,t_3,x_3,v,v,k,l_0,w,phi,'ad');
% f_4k = fK(t_4,t_3,x_3,v,v,k,l_0,w,phi,'ad');
% 
% % Affichage
% figure(1)
% subplot(3,1,1); hold on
% % regime 1 : glissement
% plot(t_01,x_01,'-','color',[0 0 1],'LineWidth',1);
% plot(t_1,x_1,'o','color',[1 0 0],'MarkerSize',3);
% % regime 2 : adherence
% plot(t_12,x_12,'-','color',[0 0 1],'LineWidth',1);
% plot(t_2,x_2,'o','color',[1 0 0],'MarkerSize',3);
% % regime 3 : glissement
% plot(t23,x23,'--','color',[0 0 0],'LineWidth',.5);
% plot(t_23,x_23,'-','color',[0 0 1],'LineWidth',1);
% plot(t_3,x_3,'o','color',[1 0 0],'MarkerSize',3)
% %% regime 4 : adherence
% %plot(t_34,x_34,'-','color',[0 0 1],'LineWidth',1);
% %plot(t_4,x_4,'o','color',[1 0 0],'MarkerSize',3);
% %
% grid('on');
% %h1 = legend('$x$','location', 'east','fontsize', 16);
% %set (h1, 'Interpreter', 'latex');
% title('Position');
% 
% subplot(3,1,2); hold on
% % regime 1 : glissement
% plot(t_01,v_01,'-','color',[.3 0 .5],'LineWidth',1);
% plot(t_1,v_1,'o','color',[1 0 0],'MarkerSize',3);
% % regime 2 : adherence
% plot(t_12,v_12,'-','color',[.3 0 .5],'LineWidth',1);
% plot(t_2,v,'o','color',[1 0 0],'MarkerSize',3);
% % regime 3 : glissement
% plot(t23,v23,'--','color',[0 0 0],'LineWidth',.5);
% plot(t_23,v_23,'-','color',[.3 0 .5],'LineWidth',1);
% plot(t_3,v_3,'o','color',[1 0 0],'MarkerSize',3);
% %% regime 4 : adherence
% %plot(t_34,v_34,'-','color',[.3 0 .5],'LineWidth',1);
% %plot(t_4,v,'o','color',[1 0 0],'MarkerSize',3);
% %
% grid('on');
% %h2 = legend('$\dot{x}$','location', 'east','fontsize', 16);
% %set (h2, 'Interpreter', 'latex');
% title('Vitesse');
% %
% subplot(3,1,3); hold on
% % regime 1 : glissement
% plot(t_01,f_01t,'-','color',[.9 .6 .7],'LineWidth',1);
% %plot([t_0,t_1],[F_c,F_c],'-','color',[.9 .6 .7],'LineWidth',1);
% plot(t_01,f_01k,'-','color',[0 .6 .2],'LineWidth',1);
% plot(t_1,f_1k,'o','color',[1 0 0],'MarkerSize',3);
% plot(t_1,f_1t,'o','color',[1 0 0],'MarkerSize',3);
% % regime 2 : adherence
% plot(t_12,f_12t,'-','color',[.9 .6 .7],'LineWidth',1);
% plot(t_12,f_12k,'-','color',[0 .6 .2],'LineWidth',1);
% plot(t_2,f_2t,'o','color',[1 0 0],'MarkerSize',1);
% plot(t_2,f_2k,'o','color',[1 0 0],'MarkerSize',3);
% % regime 3 : glissement
% plot(t23,f23k,'--','color',[0 0 0],'LineWidth',.5);
% plot(t23,f23t,'--','color',[0 0 0],'LineWidth',.5);
% %
% plot(t_23,f_23t,'-','color',[.9 .6 .7],'LineWidth',1);
% plot(t_23,f_23k,'-','color',[0 .6 .2],'LineWidth',1);
% plot(t_3,F_c,'o','color',[1 0 0],'MarkerSize',1);
% plot(t_3,f_3k,'o','color',[1 0 0],'MarkerSize',3);
% % regime 4 : adherence
% %plot(t_34,f_34t,'-','color',[.9 .6 .7],'LineWidth',1);
% %plot(t_34,f_34k,'-','color',[0 .6 .2],'LineWidth',1);
% plot([t_0,max(t23)],[F_c,F_c],':','color',[0 0 0],'LineWidth',.5);
% plot([t_0,max(t23)],[-F_c,-F_c],':','color',[0 0 0],'LineWidth',.5);
% %
% h3 = legend('$F_T$','$F_K$');%,'location', 'east','fontsize', 16);
% set (h3, 'Interpreter', 'latex');
% grid('on');
% title('Force');

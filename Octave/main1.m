close all; clear all; clc;
format short

%% Variables d'entree
[k,l_0,m,g,v,nu]=VarEntree;

%% Variables intermediaires
[F_c,w,tcF,tcK]=VarInter(k,l_0,m,g,v,nu);

%% Variables initialisation
t_0 = 0;
x_0 = l_0;
v_0 = v;
nt = 200;
%% Phase 1 : adherence ------------------------
phi = Phi(v_0,k,l_0,v,F_c);
% ---
t01 = linspace(t_0,t_0 + 4*max(tcK,tcF),nt);
x01 = xA(t01,t_0,x_0,v_0);
v01 = v*ones(size(t01));
fk01 = fK(t01,t_0,x_0,v_0,v,k,l_0,w,phi,'ad');
ft01 = fT(t01,t_0,x_0,v_0,v,k,l_0,F_c,'ad');

CostA = @(t,t_d,x_d,v_d,v,k,l_0,F_c) (abs(fT(t,t_d,x_d,v_d,v,k,l_0,F_c,'ad') - F_c).^2);
it2 = find(diff(sign(diff(CostA(t01,t_0,x_0,v_0,v,k,l_0,F_c))))==2,1);

t_1 = fminsearch(@(t) CostA(t,t_0,x_0,v_0,v,k,l_0,F_c),t01(it2+1))
t_01 = linspace(t_0,t_1,nt);
x_01 = xA(t_01,t_0,x_0,v_0);
v_01 = v*ones(size(t_01));
fk_01 = fK(t_01,t_0,x_0,v_0,v,k,l_0,w,phi,'ad');
ft_01 = fT(t_01,t_0,x_0,v_0,v,k,l_0,F_c,'ad');

%% Phase 2 : glissement ------------------------
x_1 = xA(t_1,t_0,x_0,v_0);
v_1 = v;
phi = Phi(v_1,k,l_0,v,F_c);
% ---
t12 = linspace(t_1,t_1 + 2*max(tcF,tcK),nt);
x12 = xG(t12,t_1,x_1,v_1,w,phi);
v12 = vG(t12,t_1,x_1,v_1,w,phi);
fk12 = fK(t12,t_1,x_1,v_1,v,k,l_0,w,phi,'gl');
ft12 = fT(t12,t_1,x_1,v_1,v,k,l_0,F_c,'gl');

CostG = @(t,t_d,x_d,v_d,w,phi,v) ((vG(t,t_d,x_d,v_d,w,phi) - v).^2);
it3 = find(diff(sign(diff(CostG(t12,t_1,x_1,v_1,w,phi,v))))==2,1);

t_2 = fminsearch(@(t) CostG(t,t_1,x_1,v_1,w,phi,v),t12(it3+1))
t_12 = linspace(t_1,t_2,nt);
x_12 = xG(t_12,t_1,x_1,v_1,w,phi);
v_12 = vG(t_12,t_1,x_1,v_1,w,phi);
fk_12 = fK(t_12,t_1,x_1,v_1,v,k,l_0,w,phi,'gl');
ft_12 = fT(t_12,t_1,x_1,v_1,v,k,l_0,F_c,'gl');

%% Affichage ------------------------------------
figure(1)
subplot(3,1,1); hold on
% regime 1 : adherence
plot(t_01,x_01,'-','color',[0 0 1],'LineWidth',1);
plot(t_1,x_01(200),'o','color',[1 0 0],'MarkerSize',3);
% regime 2 : glissement
plot(t_12,x_12,'-','color',[0 0 1],'LineWidth',1);
plot(t_2,x_12(200),'o','color',[1 0 0],'MarkerSize',3);
plot(t12,x12,'--','color',[0 0 0],'LineWidth',.5);
%
grid('on');
h1 = legend('$x$','location', 'east','fontsize', 16);
set (h1, 'Interpreter', 'latex');
title('Position');

subplot(3,1,2); hold on
% regime 1 : adherence
plot(t_01,v_01,'-','color',[.3 0 .5],'LineWidth',1);
plot(t_1,v_01(200),'o','color',[1 0 0],'MarkerSize',3);
% regime 2 : glissement
plot(t_12,v_12,'-','color',[.3 0 .5],'LineWidth',1);
plot(t_2,v,'o','color',[1 0 0],'MarkerSize',3);
plot(t12,v12,'--','color',[0 0 0],'LineWidth',.5);
%
grid('on');
h2 = legend('$\dot{x}$','location', 'east','fontsize', 16);
set (h2, 'Interpreter', 'latex');
title('Vitesse');
%
subplot(3,1,3); hold on
% regime 1 : adherence
plot(t_01,ft_01,'-','color',[.9 .6 .7],'LineWidth',1);
plot(t_01,fk_01,'-','color',[0 .6 .2],'LineWidth',1);
plot(t_1,fk_01(200),'o','color',[1 0 0],'MarkerSize',3);
plot(t_1,ft_01(200),'o','color',[1 0 0],'MarkerSize',3);
% regime 2 : glissement
plot(t_12,ft_12,'-','color',[.9 .6 .7],'LineWidth',1);
plot(t_12,fk_12,'-','color',[0 .6 .2],'LineWidth',1);
plot(t_2,ft_12(200),'o','color',[1 0 0],'MarkerSize',1);
plot(t_2,fk_12(200),'o','color',[1 0 0],'MarkerSize',3);
plot(t12,fk12,'--','color',[0 0 0],'LineWidth',.5);
plot(t12,ft12,'--','color',[0 0 0],'LineWidth',.5);
%
plot([t_0,max(t12)],[F_c,F_c],':','color',[0 0 0],'LineWidth',.5);
plot([t_0,max(t12)],[-F_c,-F_c],':','color',[0 0 0],'LineWidth',.5);

%
h3 = legend('$F_T$','$F_K$');%,'location', 'east','fontsize', 16);
set (h3, 'Interpreter', 'latex');
grid('on');
title('Force');

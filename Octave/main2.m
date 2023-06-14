close all; clear all; clc;
format short

%% Variables d'entree
[k,l_0,m,g,v,nu]=VarEntree;

%% Variables intermediaires
[F_c,w,tcF,tcK]=VarInter(k,l_0,m,g,v,nu);

%% Variables initialisation
t_0 = 0
x_0 = 0.4
v_0 = v;

nt=200;

%%Phase 1 : glissement ----------------------
phi = Phi(v_0,k,l_0,v,F_c);
% ---
t01=linspace(t_0,t_0 + 4*max(tcF,tcK),nt);
x01 = xG(t01,t_0,x_0,v_0,w,phi);
v01 = vG(t01,t_0,x_0,v_0,w,phi);
fk01 = fK(t01,t_0,x_0,v_0,v,k,l_0,w,phi,'gl');
ft01 = fT(t01,t_0,x_0,v_0,v,k,l_0,F_c,'gl');

CostG = @(t,t_d,x_d,v_d,w,phi,v) ((vG(t,t_d,x_d,v_d,w,phi) - v).^2);
it1=find(diff(sign(diff(CostG(t01,t_0,x_0,v_0,w,phi,v))))==2,1);

t_1 = fminsearch(@(t) CostG(t,t_0,x_0,v_0,w,phi,v),t01(it1+1))
t_01 = linspace(t_0,t_1,nt);
x_01 = xG(t_01,t_0,x_0,v_0,w,phi);
v_01 = vG(t_01,t_0,x_0,v_0,w,phi);
fk_01 = fK(t_01,t_0,x_0,v_0,v,k,l_0,w,phi,'gl');
ft_01 = fT(t_01,t_0,x_0,v_0,v,k,l_0,F_c,'gl');

%% Phase 2 : adherence ------------------------
x_1=xG(t_1,t_0,x_0,v_0,w,phi);
v_1=v;
phi = Phi(v_1,k,l_0,v,F_c);
% ---
t12 = linspace(t_1,t_1 + 4*max(tcK,tcF),nt);
x12 = xA(t12,t_1,x_1,v_1);
v12 = v*ones(size(t12));
fk12 = fK(t12,t_1,x_1,v_1,v,k,l_0,w,phi,'ad');
ft12 = fT(t12,t_1,x_1,v_1,v,k,l_0,F_c,'ad');

CostA = @(t,t_d,x_d,v_d,v,k,l_0,F_c) (abs(fT(t,t_d,x_d,v_d,v,k,l_0,F_c,'ad')) - F_c).^2;
it2 = find(diff(sign(diff(CostA(t12,t_1,x_1,v_1,v,k,l_0,F_c))))==2,1);

t_2 = fminsearch(@(t) CostA(t,t_1,x_1,v_1,v,k,l_0,F_c),t12(it2+1))
t_12 = linspace(t_1,t_2,nt);
x_12 = xA(t_12,t_1,x_1,v_1);
v_12 = v*ones(size(t_12));
fk_12 = fK(t_12,t_1,x_1,v_1,v,k,l_0,w,phi,'ad');
ft_12 = fT(t_12,t_1,x_1,v_1,v,k,l_0,F_c,'ad');

%% Phase 3 : glissement ------------------------
x_2 = xA(t_2,t_1,x_1,v_1);
v_2 = v;
phi = Phi(v_2,k,l_0,v,F_c);
% ---
t23 = linspace(t_2,t_2 + 2*max(tcF,tcK),nt);
x23 = xG(t23,t_2,x_2,v_2,w,phi);
v23 = vG(t23,t_2,x_2,v_2,w,phi);
fk23 = fK(t23,t_2,x_2,v_2,v,k,l_0,w,phi,'gl');
ft23 = fT(t23,t_2,x_2,v_2,v,k,l_0,F_c,'gl');

it3 = find(diff(sign(diff(CostG(t23,t_2,x_2,v_2,w,phi,v))))==2,1);

t_3 = fminsearch(@(t) CostG(t,t_2,x_2,v_2,w,phi,v),t23(it3+1))
t_23 = linspace(t_2,t_3,nt);
x_23 = xG(t_23,t_2,x_2,v_2,w,phi);
v_23 = vG(t_23,t_2,x_2,v_2,w,phi);
fk_23 = fK(t_23,t_2,x_2,v_2,v,k,l_0,w,phi,'gl');
ft_23 = fT(t_23,t_2,x_2,v_2,v,k,l_0,F_c,'gl');

%% Affichage ------------------------------------
figure(1)
subplot(3,1,1); hold on
% regime 1 : glissement
plot(t_01,x_01,'-','color',[0 0 1],'LineWidth',1);
plot(t_1,x_01(200),'o','color',[1 0 0],'MarkerSize',3);
% regime 2 : adherence
plot(t_12,x_12,'-','color',[0 0 1],'LineWidth',1);
plot(t_2,x_12(200),'o','color',[1 0 0],'MarkerSize',3);
% regime 3 : glissement
plot(t23,x23,'--','color',[0 0 0],'LineWidth',.5);
plot(t_23,x_23,'-','color',[0 0 1],'LineWidth',1);
plot(t_3,x_23(200),'o','color',[1 0 0],'MarkerSize',3)
%
grid('on');
h1 = legend('$x$','location', 'east','fontsize', 16);
set (h1, 'Interpreter', 'latex');
title('Position');

subplot(3,1,2); hold on
% regime 1 : glissement
plot(t_01,v_01,'-','color',[.3 0 .5],'LineWidth',1);
plot(t_1,v_01(200),'o','color',[1 0 0],'MarkerSize',3);
% regime 2 : adherence
plot(t_12,v_12,'-','color',[.3 0 .5],'LineWidth',1);
plot(t_2,v,'o','color',[1 0 0],'MarkerSize',3);
% regime 3 : glissement
plot(t23,v23,'--','color',[0 0 0],'LineWidth',.5);
plot(t_23,v_23,'-','color',[.3 0 .5],'LineWidth',1);
plot(t_3,v_23(200),'o','color',[1 0 0],'MarkerSize',3);
%
grid('on');
h2 = legend('$\dot{x}$','location', 'east','fontsize', 16);
set (h2, 'Interpreter', 'latex');
title('Vitesse');
%
subplot(3,1,3); hold on
% regime 1 : glissement
plot(t_01,ft_01,'-','color',[.9 .6 .7],'LineWidth',1);
plot(t_01,fk_01,'-','color',[0 .6 .2],'LineWidth',1);
plot(t_1,fk_01(200),'o','color',[1 0 0],'MarkerSize',3);
plot(t_1,ft_01(200),'o','color',[1 0 0],'MarkerSize',3);
% regime 2 : adherence
plot(t_12,ft_12,'-','color',[.9 .6 .7],'LineWidth',1);
plot(t_12,fk_12,'-','color',[0 .6 .2],'LineWidth',1);
plot(t_2,ft_12(200),'o','color',[1 0 0],'MarkerSize',1);
plot(t_2,fk_12(200),'o','color',[1 0 0],'MarkerSize',3);
% regime 3 : glissement
plot(t23,fk23,'--','color',[0 0 0],'LineWidth',.5);
plot(t23,ft23,'--','color',[0 0 0],'LineWidth',.5);
%
plot(t_23,ft_23,'-','color',[.9 .6 .7],'LineWidth',1);
plot(t_23,fk_23,'-','color',[0 .6 .2],'LineWidth',1);
plot(t_3,F_c,'o','color',[1 0 0],'MarkerSize',1);
plot(t_3,fk_23(200),'o','color',[1 0 0],'MarkerSize',3);
%
plot([t_0,max(t23)],[F_c,F_c],':','color',[0 0 0],'LineWidth',.5);
plot([t_0,max(t23)],[-F_c,-F_c],':','color',[0 0 0],'LineWidth',.5);
%
h3 = legend('$F_T$','$F_K$');%,'location', 'east','fontsize', 16);
set (h3, 'Interpreter', 'latex');
grid('on');
title('Force');

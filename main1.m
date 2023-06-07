close all; clear all; clc;
format long

%% Variables d'entree
[k,l_0,m,g,v,nu]=VarEntree;

%% Variables intermediaires
[F_c,omega,tcF,tcK]=VarInter(k,l_0,m,g,v,nu);
phi = Phi(v,k,l_0,v,F_c);

%% Variables initialisation
t_0 = 0     % tmps ini       [s]
x_0 = l_0;  % position ini   [m]

nt=200;

%% Phase 1 : adherence --------------------------
% Estimation numerique de t_1
Cost1=@(t) (fT(t,t_0,x_0,v,v,k,l_0,F_c,'ad')-F_c).^2;
%-------------------------------------------------
% Temps
t_1 = fminsearch(@(t) Cost1(t),t_0+tcK)
t1A =  t_0+tcF - (x_0 - l_0)/v;
ecart_t1 = t1A - t_1
t_01 = linspace(t_0,t_1,nt);
% Position
x_1 = xA(t_1,t_0,x_0,v);
x_01 = xA(t_01,t_0,x_0,v);
% Vitesse
v_01 = v*ones(size(t_01));
% Force du tapis
f_01t = fT(t_01,t_0,x_0,v,v,k,l_0,F_c,'ad');
% Force de rappel
f_01k = fK (t_01,t_0,x_0,v,v,k,l_0,omega,phi,'ad');

%% Phase 2 : glissement -------------------------
Cost2 = @(t) (d_xG(t,t_1,x_1,v,omega,phi) - v).^2;
%-------------------------------------------------
% Temps
t_2 = fminsearch(@(t) Cost2(t),t_1+tcK)
t_12 = linspace(t_1,t_2,nt);
t12 = linspace(t_1,t_1+4*max(tcF,tcK),nt);
% Position
x_2 = xG(t_2,t_1,x_1,v,v,omega,phi);
x_12 = xG(t_12,t_1,x_1,v,v,omega,phi);
x12 = xG(t12,t_1,x_1,v,v,omega,phi);
% Vitesse
v_2 = d_xG(t_2,t_1,x_1,v,omega,phi);
v_12 = d_xG(t_12,t_1,x_1,v,omega,phi);
v12 = d_xG(t12,t_1,x_1,v,omega,phi);
% Force du tapis
f_2k =  fT(t_2,t_1,x_1,v,v,k,l_0,F_c,'gl');
f_12k = fT(t_12,t_1,x_1,v,v,k,l_0,F_c,'gl');
% Force du ressort
f_2k =  fK(t_2,t_1,x_1,v,v,k,l_0,omega,phi,'gl');
f_12k = fK(t_12,t_1,x_1,v,v,k,l_0,omega,phi,'gl');

%% Phase 3 : adherence --------------------------
Cost=@(t) (fT(t,t_2,x_2,v,v,k,l_0,F_c,'ad')-F_c).^2;
%-------------------------------------------------
% Temps
t_3 = fminsearch(@(t) Cost(t),t_2+tcK);
%t3A =  t_2+tcF - (x_2 - x_0)/v;
% ecart_t3 = t3A - t_3;
t_23 = linspace(t_2,t_3,nt);
t23 = t_2+linspace(0,4*max(tcF,tcK),nt);
% Position
x_3 = xA(t_3,t_2,x_2,v);
x_23 = xA(t_23,t_2,x_2,v);
% Vitesse
v_23 = v*ones(size(t_23));
% Force tapis
f_3t = fT(t_3,t_2,x_2,v,v,k,l_0,F_c,'ad');
f_23t = fT(t_23,t_2,x_2,v,v,k,l_0,F_c,'ad');
% Force du ressort
f_3k =  fK(t_3,t_2,x_2,v,v,k,l_0,omega,phi,'ad');
f_23k = fK(t_23,t_2,x_2,v,v,k,l_0,omega,phi,'ad');


%% Affichage -------------------------------------
figure(1)

subplot(3,1,1); hold on
% regime d'adherence
plot(t_01,x_01,'-','color',[0 0 1],'LineWidth',1);
plot(t_1,x_1,'o','color',[1 0 0],'MarkerSize',3);
% regime de glissement
plot(t12,x12,'--','color',[0 0 0],'LineWidth',.5);
plot(t_12,x_12,'-','color',[0 0 1],'LineWidth',1);
plot(t_2,x_2,'o','color',[1 0 0],'MarkerSize',3);
%
grid('on');
h1 = legend('$x$','location', 'east','fontsize', 16);
set (h1, 'Interpreter', 'latex');
title('Position');

subplot(3,1,2); hold on
% regime 1 : adherence
plot([0,t_1],[v,v],'-','color',[.3 0 .5],'LineWidth',1);
plot(t_1,v,'o','color',[1 0 0],'MarkerSize',3);
% regime 2 : glissement
plot(t12,v12,'--','color',[0 0 0],'LineWidth',.5);
plot(t_12,v_12,'-','color',[.3 0 .5],'LineWidth',1);
plot(t_2,v_2,'o','color',[1 0 0],'MarkerSize',3);
%
grid('on');
h = legend('$\dot{x}_G$','location', 'east','fontsize', 16);
set (h, 'Interpreter', 'latex');
title('Vitesse');

subplot(3,1,3); hold on
% regime 1 : adherence
plot(t_01,f_01t,'-','color',[.9 .6 .7],'LineWidth',1);
plot(t_01,f_01k,'-','color',[0 .6 .2],'LineWidth',1);
plot(t_1,F_c,'o','color',[1 0 0],'MarkerSize',3);
plot(t_1,fK(t_2,t_1,x_1,v,v,k,l_0,omega,phi,'gl'),'o','color',[1 0 0],'MarkerSize',3);
% regime 2 : glissement
plot([t_2,max(t12)],[F_c,F_c],'--','color',[0 0 0],'LineWidth',.5);
plot(t12,fK(t12,t_1,x_1,v,v,k,l_0,omega,phi,'gl'),'--','color',[0 0 0],'LineWidth',.5);
plot([t_1,t_2],[F_c,F_c],'-','color',[.9 .6 .7],'LineWidth',1);
plot(t_12,f_12k,'-','color',[0 .6 .2],'LineWidth',1);
plot(t_2,F_c,'o','color',[1 0 0],'MarkerSize',1);
plot(t_2,f_2k,'o','color',[1 0 0],'MarkerSize',3);
%
grid('on');
h2 = legend('$F_T$','$F_K$','location', 'east','fontsize', 16);
set(h2, 'Interpreter', 'latex');
title('Force');
%print('figure2.pdf','-dpdf');

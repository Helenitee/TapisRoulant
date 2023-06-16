%lst1 =  List(l_0+0.1,l_0+0.5,-1,1,30,30,k,l_0,m,g,v,nu,F_c,w,tcF,tcK);
##
## x = lst1(:,1);
## y = lst1(:,2);
## z1 = lst1(:,3);
## z2 = lst1(:,6);
##
## figure(3);
## plot3(x,y,z1,'-+','color',[0 1 1]);
## xlabel('x_0'),ylabel('v_0'),zlabel('t_1'),
## grid('on');
##
##
## figure(4)
## plot3(x,y,z2,'-+','color',[1 0 1]);
## xlabel('x_0'),ylabel('v_0'),zlabel('t_2 - t_1'),
## grid('on');

## lst2 =  List(l_0-0.34,l_0-0.1004,-1,1,30,30,k,l_0,m,g,v,nu,F_c,w,tcF,tcK);
##
## x = [lst2(:,1);lst1(:,1)];
## y = [lst2(:,2);lst1(:,2)];
## z1 = [lst2(:,3);lst1(:,3)];
## z2 = [lst2(:,6);lst1(:,6)];

## figure(3);
## plot3(x,y,z1,'+','color',[0 1 1]);
## xlabel('x_0'),ylabel('v_0'),zlabel('t_1'),
## grid('on');
## title('Représentation graphique de t_1 suivant x_0 et v_0')

## figure(4)
## plot3(x,y,z2,'+','color',[1 0 1]);
## xlabel('x_0'),ylabel('v_0'),zlabel('t_2 - t_1'),
## grid('on');
## title('Représentation de t_2 - t_1 suivant x_0 et v_0')

 x=l(:,1);
 y=l(:,2);
 z=l(:,6);

 figure(7)
 surf(x,y,z,'+','color',[1 0 1]);
 xlabel('x_0'),ylabel('v_0'),zlabel('t_1'),
 grid('on');
 title('Représentation de t_1 suivant x_0 et v_0')













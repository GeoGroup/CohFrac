figure(1)
clear;clc;close all;
Globals;
rng(123456781);
rmax=12;rmin=6;
dr=0.25;   eta=0.20; n1=6;  n2=12;  ext=1.5;
rot=-1;distance=0.1;
p = particle_generate_polygon_2d(rmin,rmax,dr,eta,n1,n2,ext,rot);
%visual_poly(p);

x1=p(:,1)';
y1=p(:,2)';
x1=[x1,p(1,1)];
y1=[y1,p(1,2)];
plot(x1,y1);

figure(2)
for i = 1:length(x1)-1
    x_qian(i)=x1(i);
    y_qian(i)=y1(i);
    x_hou(i)=x1(i+1);
    y_hou(i)=y1(i+1);
end

for i = 1:length(x1)-1
    xx(i,:)=[x_qian(i),x_hou(i)];
    yy(i,:)=[y_qian(i),y_hou(i)];
end
plot(xx,yy,'k')
hold on

xx_middle = (x_qian+x_hou)/2;
yy_middle = (y_qian+y_hou)/2;
scatter(xx_middle,yy_middle,'r')

figure(3)
rand_max = 0.7;
rand_min = 0.3;
for i =1 :length(xx_middle)
   set1_len(i) = sqrt((x_qian(i)-x_hou(i))^2+(y_qian(i)-y_hou(i))^2);
end

amplen=0.05*set1_len*rand;


vert_x = x_qian-x_hou;
vert_y = y_qian-y_hou;
vert = [vert_x,vert_y];
vert_norm_x = vert_x./sqrt(vert_x.^2+vert_y.^2);
vert_norm_y = vert_y./sqrt(vert_x.^2+vert_y.^2);
vert_norm = [vert_norm_x,vert_norm_y]; %正交化
rot_1=deg2rad(90);
vert_x_r1 = vert_norm_x.*cos(rot_1)-vert_norm_y.*sin(rot_1);
vert_y_r1 = vert_norm_x.*sin(rot_1)+vert_norm_y.*cos(rot_1);
vert_r1 = [vert_x_r1,vert_y_r1];

x_amp = xx_middle+vert_x_r1.*amplen;
y_amp = yy_middle+vert_y_r1.*amplen;

for i =1 :length(xx_middle)
    xxx(i,:) = [x_qian(i),x_amp(i),x_hou(i)];
    yyy(i,:) = [y_qian(i),y_amp(i),y_hou(i)];
    plot(xxx(i,:),yyy(i,:),'r');
    hold on
end

figure(4)
L=4;     %%%%% L表示分形迭代次数
d=0.5;
for i = 1 :length(xx_middle)
    [xpp1(i,:),ypp1(i,:),m]=iteration_of_fractal(xxx(i,:),yyy(i,:),d,L);%%%%%%%%%%%%%%%%%%%
    for j = 1:size(xpp1(i,:))
        plot(xpp1(i,:),ypp1(i,:),'k-');
        hold on
    end
end
% figure(5)
% 
% k=1;
% for i = 2:length(xx_middle)-1
%    ab=[xpp1(i,:)-xpp1(i-1,:),ypp1(i,:)-ypp1(i-1,:)];
%    bc=[xpp1(i+1,:)-xpp1(i,:),ypp1(i+1,:)-ypp1(i,:)];
%    sigma(i) = acos(dot(ab,bc)/(norm(ab)*norm(bc)));
%    sigma(i) = sigma(i)*180/pi;
%    if sigma(i) < 100
%        mm(k,:)=xpp1(i,:);
%        nn(k,:)=ypp1(i,:);
%        k = k+1;
%    end
% end
% 
% plot(mm,nn)





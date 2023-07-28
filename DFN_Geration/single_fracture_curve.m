
clear;clc;close all;
Globals;
rng(1234567890);

set1 = [0.7 0.5-0.1*sqrt(3) 0.4 0.5+0.2*sqrt(3)];

p1=[0.5,0.48];p2=[0.5,0.52];
plot([p1(1),p2(1)],[p1(2),p2(2)],'LineWidth',2,'color','b');hold on
 %初始图像，线段角度固定
figure(1)
Draw('lin',set1);
axis([0 1 0 1])
%% 
%线段在指定范围内随机旋转
figure(2)

x1 = set1(:,1);
y1 = set1(:,2);
x2 = set1(:,3);
y2 = set1(:,4);
[set1_row,set1_column] =size(x1);
x_middle = (x1+x2)/2;
y_middle = (y1+y2)/2;
scatter(x_middle,y_middle,'r')
hold on
theta_max=deg2rad(0);
theta = -theta_max+2*theta_max*rand([set1_row,1]);
x1_r = (x1-x_middle).*cos(theta)-(y1-y_middle).*sin(theta)+x_middle;
y1_r = (x1-x_middle).*sin(theta)+(y1-y_middle).*cos(theta)+y_middle;
x2_r = (x2-x_middle).*cos(theta)-(y2-y_middle).*sin(theta)+x_middle;
y2_r = (x2-x_middle).*sin(theta)+(y2-y_middle).*cos(theta)+y_middle;
Draw('lin',[x1_r,y1_r,x2_r,y2_r]);
hold on
%% 
%将直线改为有振幅的折线
figure(3)
%在直线上随机取两点，沿垂直方向随机延伸一定距离
%x3=(x1_r+x_middle)*rand
rand_max = 0.7;
rand_min = 0.3;
amp_fac_1 = rand_min+(rand_max-rand_min)*rand([set1_row,1]);  %振幅因子
amp_fac_2 = rand_min+(rand_max-rand_min)*rand([set1_row,1]);

x3 = x1_r.*amp_fac_1+x_middle.*(ones(set1_row,1)-amp_fac_1);
y3 = y1_r.*amp_fac_1+y_middle.*(ones(set1_row,1)-amp_fac_1);
x4 = x2_r.*amp_fac_2+x_middle.*(ones(set1_row,1)-amp_fac_2);
y4 = y2_r.*amp_fac_2+y_middle.*(ones(set1_row,1)-amp_fac_2);

scatter(x3,y3,'b');
hold on
scatter(x4,y4,'b');
hold on
Draw('lin',[x1_r,y1_r,x2_r,y2_r]);
hold on

%% 
%选中的两个点向外延伸
figure(4)
vert_x = x1_r-x2_r;
vert_y = y1_r-y2_r;
vert = [vert_x,vert_y];
vert_norm_x = vert_x./sqrt(vert_x.^2+vert_y.^2);
vert_norm_y = vert_y./sqrt(vert_x.^2+vert_y.^2);
vert_norm = [vert_norm_x,vert_norm_y]; %正交化

rot_1=deg2rad(90);
rot_2=deg2rad(-90);
vert_x_r1 = vert_norm_x.*cos(rot_1)-vert_norm_y.*sin(rot_1);
vert_y_r1 = vert_norm_x.*sin(rot_1)+vert_norm_y.*cos(rot_1);
vert_r1 = [vert_x_r1,vert_y_r1];
vert_x_r2 = vert_norm_x.*cos(rot_2)-vert_norm_y.*sin(rot_2);
vert_y_r2 = vert_norm_x.*sin(rot_2)+vert_norm_y.*cos(rot_2);
vert_r2 = [vert_x_r2,vert_y_r2];

for i =1 :set1_row
   set1_len(i) = sqrt((x1(i)-x2(i))^2+(y1(i)-y2(i))^2);
end
set1_len = set1_len';
amplen_1=0.05*set1_len.*rand([set1_row,1]);
amplen_2=0.05*set1_len.*rand([set1_row,1]);
x_amp_1 = x3+vert_x_r1.*amplen_1;
y_amp_1 = y3+vert_y_r1.*amplen_1;
x_amp_2 = x4+vert_x_r2.*amplen_2;
y_amp_2 = y4+vert_y_r2.*amplen_2;
% Draw('lin',[x_amp_1,y_amp_1,x_amp_2,y_amp_2]);
% hold on
% Draw('lin',[x1_r,y1_r,x_amp_1,y_amp_1]);
% hold on
% Draw('lin',[x_amp_2,y_amp_2,x2_r,y2_r]); 
xx = [x1_r,x_amp_1,x_amp_2,x2_r];
yy = [y1_r,y_amp_1,y_amp_2,y2_r];
for i = 1:set1_row
    plot(xx(i,:),yy(i,:),'LineWidth',2,'color','r');
    hold on
end
axis([0 1 0 1])
hold on

%% 
figure(5)
L=4;     %%%%% L表示分形迭代次数
d=0.45;

% [xpp,ypp,m]=iteration_of_fractal(xx(1,:),yy(1,:),d,L);%%%%%%%%%%%%%%%%%%%
% for j = 1:size(xpp)
%     plot(xpp(j,:),ypp(j,:),'k-');
%     hold on
% end

fid1=fopen('single.dxf','w');
write_dxf_head(fid1);

for i =1:set1_row
    [xpp,ypp,m_1]=iteration_of_fractal(xx(i,:),yy(i,:),d,L);
    plot(xpp,ypp,'LineWidth',2,'color','r');
    hold on
    line1(:,1)=xpp(1:m_1-2)';
    line1(:,2)=ypp(1:m_1-2)';
    line1(:,3)=0;
    line1(:,4)=xpp(2:m_1-1)';
    line1(:,5)=ypp(2:m_1-1)';
    line1(:,6)=0; 
    write_dxf_line(fid1,line1);
end
axis([0 1 0 1])
hold on
write_dxf_end(fid1);


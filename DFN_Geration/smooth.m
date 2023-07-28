%%%%%%%%%%%%%Modified DFN%%%%%%%%%%%%
%% 
clear;clc;close all;
Globals;
rng(1234567890);
line_num = 100;
line_num_2 = 50;
set1 = Field(DFN('dim',2,'n',line_num,'dir',0,'ddir',-1e9,'minl',0.2,...
            'mu',0.2,'maxl',0.2,'dsep',0.001,'bbx',[0,0,1,1]),'Line');
set2 = Field(DFN('dim',2,'n',line_num_2,'dir',90,'ddir',-1e9,'minl',0.2,...
            'mu',0.2,'maxl',0.2,'dsep',0.001,'bbx',[0,0,1,1]),'Line');

 %初始图像，线段角度固定
figure(1)
Draw('lin',[set1;set2]);
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

m1 = set2(:,1);
n1 = set2(:,2);
m2 = set2(:,3);
n2 = set2(:,4);
[set2_row,set2_column] =size(m1);
m_middle = (m1+m2)/2;
n_middle = (n1+n2)/2;
scatter(m_middle,n_middle,'r');
hold on
theta_max_2=deg2rad(0);
theta_2 = -theta_max_2+2*theta_max_2*rand([set2_row,1]);
m1_r = (m1-m_middle).*cos(theta_2)-(n1-n_middle).*sin(theta_2)+m_middle;
n1_r = (m1-m_middle).*sin(theta_2)+(n1-n_middle).*cos(theta_2)+n_middle;
m2_r = (m2-m_middle).*cos(theta_2)-(n2-n_middle).*sin(theta_2)+m_middle;
n2_r = (m2-m_middle).*sin(theta_2)+(n2-n_middle).*cos(theta_2)+n_middle;

Draw('lin',[m1_r,n1_r,m2_r,n2_r]);
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

amp_fac_3 = rand_min+(rand_max-rand_min)*rand([set2_row,1]);  %振幅因子
amp_fac_4 = rand_min+(rand_max-rand_min)*rand([set2_row,1]);

m3 = m1_r.*amp_fac_3+m_middle.*(ones(set2_row,1)-amp_fac_3);
n3 = n1_r.*amp_fac_3+n_middle.*(ones(set2_row,1)-amp_fac_3);
m4 = m2_r.*amp_fac_4+m_middle.*(ones(set2_row,1)-amp_fac_4);
n4 = n2_r.*amp_fac_4+n_middle.*(ones(set2_row,1)-amp_fac_4);

scatter(m3,n3,'b');
hold on
scatter(m4,n4,'b');
hold on
Draw('lin',[m1_r,n1_r,m2_r,n2_r]);
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
amplen_1=0.4*set1_len.*rand([set1_row,1]);
amplen_2=0.4*set1_len.*rand([set1_row,1]);
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

vert_m = m1_r-m2_r;
vert_n = n1_r-n2_r;
vert_2 = [vert_m,vert_n];
vert_norm_m = vert_m./sqrt(vert_m.^2+vert_n.^2);
vert_norm_n = vert_n./sqrt(vert_m.^2+vert_n.^2);
vert_norm_2 = [vert_norm_m,vert_norm_n]; %正交化

rot_3=deg2rad(90);
rot_4=deg2rad(-90);
vert_m_r1 = vert_norm_m.*cos(rot_3)-vert_norm_n.*sin(rot_3);
vert_n_r1 = vert_norm_m.*sin(rot_3)+vert_norm_n.*cos(rot_3);
vert_r1_2 = [vert_m_r1,vert_n_r1];
vert_m_r2 = vert_norm_m.*cos(rot_4)-vert_norm_n.*sin(rot_4);
vert_n_r2 = vert_norm_m.*sin(rot_4)+vert_norm_n.*cos(rot_4);
vert_r2_2 = [vert_m_r2,vert_n_r2];

for i =1 :set2_row
   set2_len(i) = sqrt((m1(i)-m2(i))^2+(n1(i)-n2(i))^2);
end
set2_len = set2_len';
amplen_3=0.4*set2_len.*rand([set2_row,1]);
amplen_4=0.4*set2_len.*rand([set2_row,1]);
m_amp_1 = m3+vert_m_r1.*amplen_3;
n_amp_1 = n3+vert_n_r1.*amplen_3;
m_amp_2 = m4+vert_m_r2.*amplen_4;
n_amp_2 = n4+vert_n_r2.*amplen_4;
% Draw('lin',[x_amp_1,y_amp_1,x_amp_2,y_amp_2]);
% hold on
% Draw('lin',[x1_r,y1_r,x_amp_1,y_amp_1]);
% hold on
% Draw('lin',[x_amp_2,y_amp_2,x2_r,y2_r]); 
mm = [m1_r,m_amp_1,m_amp_2,m2_r];
nn = [n1_r,n_amp_1,n_amp_2,n2_r];
for i = 1:set2_row
    plot(mm(i,:),nn(i,:),'LineWidth',2,'color','r');
    hold on
end
axis([0 1 0 1])
%% 
figure(5)
L=4;     %%%%% L表示分形迭代次数
d=0.4;

% [xpp,ypp,m]=iteration_of_fractal(xx(1,:),yy(1,:),d,L);%%%%%%%%%%%%%%%%%%%
% for j = 1:size(xpp)
%     plot(xpp(j,:),ypp(j,:),'k-');
%     hold on
% end

fid1=fopen('smooth.dxf','w');
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

for i =1:set2_row
    [mpp,npp,m_2]=iteration_of_fractal(mm(i,:),nn(i,:),d,L);%%%%%%%%%%%%%%%%%%%
    for j = 1:size(mpp,1)
        plot(mpp,npp,'LineWidth',2,'color','r');
        hold on
        line2(:,1)=mpp(1:m_2-2)';
        line2(:,2)=npp(1:m_2-2)';
        line2(:,3)=0;
        line2(:,4)=mpp(2:m_2-1)';
        line2(:,5)=npp(2:m_2-1)';
        line2(:,6)=0; 
        write_dxf_line(fid1,line2);
    end
end
write_dxf_end(fid1);
%% 
% figure(6)
% d=[0.1 0.3 0.5 0.7];
% L=4;     %%%%% L表示分形迭代次数
% for q=1:4 
%     for i = 1:line_num
%         [xpp,ypp,m_1]=iteration_of_fractal(xx(i,:),yy(i,:),d(q),L);%%%%%%%%%%%%%%%%%%%
%         subplot(2,2,q);
%         for j = 1:size(xpp)
%             plot(xpp(j,:),ypp(j,:),'LineWidth',2,'color','r');
%             hold on
%         end
%     end
% end
% hold on
% 
% for q=1:4 
%     for i = 1:set2_row
%         [mpp,npp,m_2]=iteration_of_fractal(mm(i,:),nn(i,:),d(q),L);%%%%%%%%%%%%%%%%%%%
%         subplot(2,2,q);
%         for j = 1:size(mpp)
%             plot(mpp(j,:),npp(j,:),'LineWidth',2,'color','r');
%             hold on
%         end
%     end
% end
% 


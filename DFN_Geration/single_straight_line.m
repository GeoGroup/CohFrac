
clear;clc;close all;
Globals;
rng(1234567890);

set1 = [0.714 0.5-0.114*sqrt(3) 0.386 0.5+0.214*sqrt(3)];

p1=[0.5,0.48];p2=[0.5,0.52];
plot([p1(1),p2(1)],[p1(2),p2(2)],'LineWidth',2,'color','b');hold on
 %初始图像，线段角度固定
figure(1)
Draw('lin',set1);
axis([0 1 0 1])


x1 = set1(:,1);
y1 = set1(:,2);
x2 = set1(:,3);
y2 = set1(:,4);
[set1_row,set1_column] =size(x1);

fid1=fopen('single_straight.dxf','w');
write_dxf_head(fid1);

line1(:,1)=0.714;
line1(:,2)=0.5-0.114*sqrt(3);
line1(:,3)=0;
line1(:,4)=0.386;
line1(:,5)=0.5+0.214*sqrt(3);
line1(:,6)=0; 

write_dxf_line(fid1,line1);

write_dxf_end(fid1);
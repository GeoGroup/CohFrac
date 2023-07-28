% %% 
clear;clc;close all;
% figure(1)
% r=10;
% theta=sort(rand(1,6)*2*pi);
% x1=r*cos(theta) ;
% x2=r*cos(theta(1));
% x=[x1 x2];
% y1=r*sin(theta) ;
% y2=r*sin(theta(1));
% y=[y1 y2];
figure(5)
w=[0 2 0.9 -1.2 0 ];
z=[-0.4 0 1.2 1 -0.4];
plot(w,z)

figure(1)
x1=[0 1 2];
y1=[-0.4 -0.5 0];
x2=[2 1.75 0.9];
y2=[0 0.5 1.2];
x3=[0.9 0 -1.2];
y3=[1.2 1.5 1];
x4=[-1.2 -1 0];
y4=[1 0.25 -0.4];

plot(x1,y1);
hold on
plot(x2,y2);
hold on
plot(x3,y3);
hold on
plot(x4,y4);
%% 
figure(2)
L=4;     %%%%% L表示分形迭代次数
d=0.4;

[xpp1,ypp1,m]=iteration_of_fractal(x1,y1,d,L);%%%%%%%%%%%%%%%%%%%
for j = 1:size(xpp1)
    plot(xpp1,ypp1,'k-');
    hold on
end

[xpp2,ypp2,m]=iteration_of_fractal(x2,y2,d,L);%%%%%%%%%%%%%%%%%%%
for j = 1:size(xpp2)
    plot(xpp2,ypp2,'k-');
    hold on
end

[xpp3,ypp3,m]=iteration_of_fractal(x3,y3,d,L);%%%%%%%%%%%%%%%%%%%
for j = 1:size(xpp3)
    plot(xpp3,ypp3,'k-');
    hold on
end

[xpp4,ypp4,m]=iteration_of_fractal(x4,y4,d,L);%%%%%%%%%%%%%%%%%%%
for j = 1:size(xpp4)
    plot(xpp4,ypp4,'k-');
    hold on
end

figure(3)
xx=[xpp1,xpp2,xpp3,xpp4,xpp1(1)];
yy=[ypp1,ypp2,ypp3,ypp4,ypp1(1)];
plot(xx,yy)
k=1
for i = 2:length(xx)-1
   ab=[xx(i)-xx(i-1),yy(i)-yy(i-1)];
   bc=[xx(i+1)-xx(i),yy(i+1)-yy(i)];
   sigma(i) = acos(dot(ab,bc)/(norm(ab)*norm(bc)));
   sigma(i) = sigma(i)*180/pi;
   if sigma(i) < 10
       mm(k)=xx(i);
       nn(k)=yy(i);
       k = k+1;
   end
end

mm=[mm,xpp1(1),xpp1(2)];
nn=[nn,ypp1(1),ypp1(2)];
plot(mm,nn)







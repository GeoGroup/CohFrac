clear;clc;close all;
figure(1)
%%裂缝长度%%
I1=xlsread('D:\MyPaper\A novel IFS-based fractal discrete fracture network for fluid-driven fracture behavior of rock mass\DFNpaper\IFS Revised manuscript\Revision3.0\数据\figure\kgd_LENGTH.xlsx');
I2=xlsread('D:\MyPaper\A novel IFS-based fractal discrete fracture network for fluid-driven fracture behavior of rock mass\DFNpaper\IFS Revised manuscript\Revision3.0\数据\figure\kgd_WIDTH.xlsx');
I3=xlsread('D:\MyPaper\A novel IFS-based fractal discrete fracture network for fluid-driven fracture behavior of rock mass\DFNpaper\IFS Revised manuscript\Revision3.0\数据\figure\kgd_POR.xlsx');

x1=I1(:,1)-1;
y1=I1(:,2)-0.1;
x2=I2(:,1)-1;
y2=I2(:,2)-0.002;
x3=I3(:,1)-1;
y3=I3(:,2);

% Timon论文的参数
Young_modulus=17e9;
Poisson_ratio=0.2;
Viscosity=0.1;
Injection_rate=0.001;

%自己论文的参数
% Young_modulus=15e9;
% Poisson_ratio=0.25;
% Viscosity=0.1;
% Injection_rate=0.001;

Strain_modulus=Young_modulus/(2*(1+Poisson_ratio));

t=linspace(0,20,200);
length=0.68*t.^(2/3).*((Strain_modulus*(Injection_rate^3))/(Viscosity*(1-Poisson_ratio)))^(1/6);
plot(t,length,'b','LineWidth',2);hold on
plot(x1,y1,'r','LineWidth',2);hold on

legend('analytical','FEM','FontSize',12,'FontName','Times New Roman','Location','NorthWest')

xlabel('t (s)','FontSize',12,'FontName','Times New Roman')
ylabel('{\it{l}}(t) (m)','FontSize',12,'FontName','Times New Roman')

axis([0 20 0 14]);
lgd=legend;
% set(lgd,'Box','off')
lgd.NumColumns = 1;
grid on

set(gca,'tickdir','out','FontName','Times New Roman','FontSize',12,'linewidth',1);
box off  
ax2 = axes('Position',get(gca,'Position'),...  
           'XAxisLocation','top',...  
           'YAxisLocation','right',...  
           'Color','none',...  
           'XColor','k','YColor','k');  
set(ax2,'YTick', [],'linewidth',1);  
set(ax2,'XTick', [],'linewidth',1);  
box on

figure(2)
%%裂缝宽度%%
width=1.87*t.^(1/3).*(Viscosity*(1-Poisson_ratio)*(Injection_rate^3)/Strain_modulus).^(1/6);
width=width*1000;
y2=y2*1000;
plot(t,width,'b','LineWidth',2);hold on
plot(x2,y2,'r','LineWidth',2);hold on

legend('analytical','FEM','FontSize',12,'FontName','Times New Roman','Location','NorthWest')

xlabel('t (s)','FontSize',12,'FontName','Times New Roman')
ylabel('{\it{w}}(0,t) (mm)','FontSize',12,'FontName','Times New Roman')

axis([0 20 0 3]);
lgd=legend;
% set(lgd,'Box','off')
lgd.NumColumns = 1;
grid on

set(gca,'tickdir','out','FontName','Times New Roman','FontSize',12,'linewidth',1);
box off  
ax2 = axes('Position',get(gca,'Position'),...  
           'XAxisLocation','top',...  
           'YAxisLocation','right',...  
           'Color','none',...  
           'XColor','k','YColor','k');  
set(ax2,'YTick', [],'linewidth',1);  
set(ax2,'XTick', [],'linewidth',1);  
box on

figure(3)
%%压力%%
pressure=1.38*((Strain_modulus^3*Injection_rate*Viscosity/((1-Poisson_ratio)^3)./length./length)).^(1/4);

pressure=pressure/1e6;
y3=y3/1e6;

plot(t,pressure,'b','LineWidth',2);hold on
plot(x3,y3,'r','LineWidth',2);hold on

legend('analytical','FEM','FontSize',12,'FontName','Times New Roman','Location','Best')

xlabel('t (s)','FontSize',12,'FontName','Times New Roman')
ylabel('{\it{p}}_{\it{f}}(0,t) (MPa)','FontSize',12,'FontName','Times New Roman')

axis([0 20 0 8]);
lgd=legend;
% set(lgd,'Box','off')
lgd.NumColumns = 1;
grid on

set(gca,'tickdir','out','FontName','Times New Roman','FontSize',12,'linewidth',1);
box off  
ax2 = axes('Position',get(gca,'Position'),...  
           'XAxisLocation','top',...  
           'YAxisLocation','right',...  
           'Color','none',...  
           'XColor','k','YColor','k');  
set(ax2,'YTick', [],'linewidth',1);  
set(ax2,'XTick', [],'linewidth',1);  
box on
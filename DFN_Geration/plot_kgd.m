clear;clc;close all;
figure(1)
%%裂缝长度%%

Young_modulus=15e9;
Poisson_ratio=0.25;
Viscosity=0.1;
Injection_rate=0.001;
% 
% Young_modulus=17e9;
% Poisson_ratio=0.2;
% Viscosity=0.1;
% Injection_rate=0.001;
Strain_modulus=Young_modulus/(2*(1+Poisson_ratio));

t=linspace(0,20,200);
length=0.68*t.^(2/3).*((Strain_modulus*(Injection_rate^3))/(Viscosity*(1-Poisson_ratio)))^(1/6);
plot(t,length);hold on

figure(2)
%%裂缝宽度%%
width=1.87*t.^(1/3).*(Viscosity*(1-Poisson_ratio)*(Injection_rate^3)/Strain_modulus).^(1/6);
plot(t,width);hold on

figure(3)
%%压力%%
pressure=1.38*((Strain_modulus^3*Injection_rate*Viscosity/((1-Poisson_ratio)^3)./length./length)).^(1/4);
plot(t,pressure);hold on

axis([0 20 0 1e7])

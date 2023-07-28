clear;clc;
x=[0 2 3 4];y=[1 0 1.5 1];
d=0.6;

L=4;   
[xpp,ypp,m]=iteration_of_fractal(x,y,d,L);%%%%%%%%%%%%%%%%%%%
plot(xpp,ypp,'k-');

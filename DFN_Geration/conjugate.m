clear;clc;close all;
Globals;
rng(1234567890);
line_num = 150;
line_num_2 = 50;
set1 = Field(DFN2('dim',2,'n',line_num,'dir',60,'ddir',-1e9,'minl',0.1,...
            'mu',0.1,'maxl',0.2,'bbx',[0,0,1,1]),'Line');
set2 = Field(DFN2('dim',2,'n',line_num_2,'dir',120,'ddir',-1e9,'minl',0.4,...
            'mu',0.2,'maxl',0.6,'bbx',[0,0,1,1]),'Line');
%初始图像，线段角度固定
figure(1)
Draw('lin',[set1;set2]);
%% 
% figure(2)


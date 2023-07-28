clc;clear;

data = xlsread('D:\abaqus2021\temp\Job-1-COH-Infos(PART-1-ALLELEMEDGES-1).csv');  
data2 = xlsread('D:\abaqus2021\temp\Job-1-COH-Output(PART-1-ALLELEMEDGES-1).csv');
data3 = importdata('C:\Users\小黄鸭\Desktop\DFN\suodian\data.txt');
crack_elem = data2(end,:);
crack_elem = crack_elem(:,6:end);
real_crack = [];
k=1;

for i = 1:length(data)
    for j = 1:length(crack_elem)
       if data(i,1) == crack_elem(1,j)
           real_crack(k,:) = data(i,:);
           k = k+1;
       end
    end
end

x1 = real_crack(:,2);
x2 = real_crack(:,4); 
y1 = real_crack(:,3); 
y2 = real_crack(:,5);
x = [x1,x2];
y = [y1,y2];
x_len = size(x1,1);
for i = 1:x_len
    plot(x(i,:),y(i,:),'LineWidth',2,'color','r')
    hold on
end


xx1 = data3(:,1);
yy1 = data3(:,2);
xx2 = data3(:,3);
yy2 = data3(:,4);
xx = [xx1,xx2];
yy = [yy1,yy2];
xx_len = size(xx1,1);
for i = 1:xx_len
    plot(xx(i,:),yy(i,:),'k')
    hold on
end



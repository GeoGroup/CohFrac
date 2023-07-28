clear;clc;
x=[0 1 2 3 4];y=[1 0.5 0 1.5 1];
d=[0.1 0.3 0.5 0.7];
L=4;     %%%%% L表示分形迭代次数
for q=1:4 
% q=5;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [xpp,ypp,m]=iteration_of_fractal(x,y,d(q),L);%%%%%%%%%%%%%%%%%%%
%     line(q,:,1)=xpp(1:m-2)';
%     line(q,:,2)=ypp(1:m-2)';
%     line(q,:,3)=0;
%     line(q,:,4)=xpp(2:m-1)';
%     line(q,:,5)=ypp(2:m-1)';
%     line(q,:,6)=0; 
%     
     subplot(2,2,q);
%     plot(x,y,'ro');
%     hold on
    plot(xpp,ypp,'k-');
    xlabel('x');
    ylabel('y');
    
end

%画边框
% xf=[0 4 4 4 0 0 0]';
% yf=[0 0 2 4 4 2 0]';
% frame(:,1)=xf(1:4);
% frame(:,2)=yf(1:4);
% frame(:,3)=0;
% frame(:,4)=xf(2:5);
% frame(:,5)=yf(2:5);
% frame(:,6)=0;
% 
fid1=fopen('linedn1.dxf','w');
line1(:,:)=line(1,:,:);
write_dxf_head(fid1);
write_dxf_line(fid1,line1);
write_dxf_end(fid1);
fid3=fopen('linedn3.dxf','w');
line3(:,:)=line(2,:,:);
write_dxf_head(fid3);
write_dxf_line(fid3,line3);
write_dxf_end(fid3);
fid5=fopen('linedn5.dxf','w');
line5(:,:)=line(3,:,:);
write_dxf_head(fid5);
write_dxf_line(fid5,line5);
write_dxf_end(fid5);
fid7=fopen('linedn7.dxf','w');
line7(:,:)=line(4,:,:);
write_dxf_head(fid7);
write_dxf_line(fid7,line7);
write_dxf_end(fid7);
% fid9=fopen('line9.dxf','w');
% line9(:,:)=line(5,:,:);
% write_dxf_head(fid9);
% write_dxf_line(fid9,line9);
% write_dxf_end(fid9);

%%%%----------Write to Gmsh----------%%%%
% fid=fopen('fractal_8_gmsh.geo','wt+');%%%%%%%%%%%%%%%%%%%%%%5
% write_geo(fid,xpp,ypp,xf,yf);
% fclose(fid1);

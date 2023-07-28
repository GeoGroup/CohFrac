function write_geo(fid,xpp,ypp,xf,yf)
% 已知结点和外边框坐标，写成GMSH可读的.geo文件

fprintf(fid,'lc1=0.1;\n'); % lc为结点网格尺寸，lc越大网格越稀疏
fprintf(fid,'lc2=0.3;\n');

% 第1个面，up面
% 转换平面方向，位于xz平面，y=0
for i=1:length(xpp)
    fprintf(fid,'%s\n',['Point(',num2str(i),')={',num2str(xpp(i)),',0,',num2str(ypp(i)),',lc1};']);
end
fprintf(fid,'%s\n',['Point(',num2str(length(xpp)+1),')={',num2str(xf(4)),',0,',num2str(yf(4)),',lc2};']);
fprintf(fid,'%s\n',['Point(',num2str(length(xpp)+2),')={',num2str(xf(5)),',0,',num2str(yf(5)),',lc2};']);

for i=1:length(xpp)+1
        fprintf(fid,'%s\n',['Line(',num2str(i),')={',num2str(i),',',num2str(i+1),'};']);
end
fprintf(fid,'%s\n',['Line(',num2str(length(xpp)+2),')={',num2str(length(xpp)+2),',1 };']);
fprintf(fid,'%s',['Line Loop(1)={']);
for i=1:length(xpp)+1
    fprintf(fid,'%s',[num2str(i),',']);
end
fprintf(fid,'%s\n',[num2str(length(xpp)+2),'};']);
fprintf(fid,'%s\n',['Plane Surface(1)={1};']);
fprintf(fid,'%s\n',['Physical Surface("up")={1};']);
% 第2个面，down面，位于xz平面，y=0
fprintf(fid,'%s\n',['Point(',num2str(length(xpp)+3),')={',num2str(xf(2)),',0,',num2str(yf(2)),',lc2};']);
fprintf(fid,'%s\n',['Point(',num2str(length(xpp)+4),')={',num2str(xf(1)),',0,',num2str(yf(1)),',lc2};']);
fprintf(fid,'%s\n',['Line(',num2str(length(xpp)+3),')={',num2str(length(xpp)),',',num2str(length(xpp)+3),' };']);
fprintf(fid,'%s\n',['Line(',num2str(length(xpp)+4),')={',num2str(length(xpp)+3),',',num2str(length(xpp)+4),' };']);
fprintf(fid,'%s\n',['Line(',num2str(length(xpp)+5),')={',num2str(length(xpp)+4),',1 };']);
fprintf(fid,'%s',['Line Loop(2)={']);
for i=1:length(xpp)-1
    fprintf(fid,'%s',[num2str(i),',']);
end
fprintf(fid,'%s\n',[num2str(length(xpp)+3),',',num2str(length(xpp)+4),',',num2str(length(xpp)+5),'};']);
fprintf(fid,'%s\n',['Plane Surface(2)={2};']);
fprintf(fid,'%s\n',['Physical Surface("down")={2};']);
end
clear;clc;close all;
Globals;
rng(1234567);
line_num = 200;
line_num_2 = 50;
set1 = Field(DFN('dim',2,'n',line_num,'dir',60,'ddir',-1e9,'minl',0.1,...
            'mu',0.1,'maxl',0.2,'bbx',[0,0,1,1]),'Line');
set2 = Field(DFN('dim',2,'n',line_num_2,'dir',120,'ddir',-1e9,'minl',0.4,...
            'mu',0.2,'maxl',0.6,'bbx',[0,0,1,1]),'Line');
%初始图像，线段角度固定
figure(1)
Draw('lin',[set1;set2]);

figure(2)
%控制裂缝间最小距离
nset1=set1(1,:);
for i = 2:size(set1,1)
   for j = 1:size(nset1,1)
       dis(j) = line_distance(set1(i,:),nset1(j,:));
   end
   min_dis = min(dis);
   if min_dis > 0.015
      nset1 = [nset1;set1(i,:)];
   end
end

nset2=set2(1,:);
for i = 2:size(set2,1)
   for j = 1:size(nset2,1)
       dis2(j) = line_distance(set2(i,:),nset2(j,:));
   end
   min_dis2 = min(dis2);
   if min_dis2 > 0.02
      nset2 = [nset2;set2(i,:)];
   end
end
Draw('lin',[nset1;nset2])
%% 
figure(3)
Draw('lin',[nset1;nset2]); hold on
trans=[];
for i = 1:size(nset1,1)
    for j = 1:size(nset2,1)
        XY1 = [nset2(j,1);nset2(j,2)];
        XY2 = [nset2(j,3);nset2(j,4)];
        out=linexlines2D(XY1,XY2,[nset1(i,1),nset1(i,2)],[nset1(i,3),nset1(i,4)]);
    %     plot(out(1,:),out(2,:),'o','MarkerFaceColor','r');hold on
        out2 = out(~isnan(out));
        [o,p]=size(out2);
        out2 = reshape(out2,2,o*p/2);
        shou = [nset1(i,1),nset1(i,2)];
        wei = [nset1(i,3),nset1(i,4)];
        if isempty(out2) == 0
            for j = 1:size(out2,2)
                trans(j,1)=out2(1,j);
                trans(j,2)=out2(2,j);
            end
            plot(trans(:,1),trans(:,2),'o','MarkerFaceColor','b'); hold on
            for k = 1:size(trans,1)
               shou_dis = norm(shou-[trans(k,1),trans(k,2)]);
               if shou_dis < 0.04
                   nset1(i,1) = trans(k,1);
                   nset1(i,2) = trans(k,2);
               end
               wei_dis = norm(wei-[trans(k,1),trans(k,2)]);
               if wei_dis < 0.04
                   nset1(i,3) = trans(k,1);
                   nset1(i,4) = trans(k,2);
               end
            end
        end
    end
end    

%%
figure(4)
Draw('lin',[nset1;nset2]);

fid1=fopen('pingban.dxf','w');
write_dxf_head(fid1);
for i =1:size(nset1,1)
    plot([nset1(i,1),nset1(i,3)],[nset1(i,2),nset1(i,4)],'LineWidth',2,'color','r');
    hold on
    line1(:,1)=nset1(i,1);
    line1(:,2)=nset1(i,2);
    line1(:,3)=0;
    line1(:,4)=nset1(i,3);
    line1(:,5)=nset1(i,4);
    line1(:,6)=0; 
    write_dxf_line(fid1,line1);
end
axis([0 1 0 1])
hold on

for i =1:size(nset2,1)
    plot([nset2(i,1),nset2(i,3)],[nset2(i,2),nset2(i,4)],'LineWidth',2,'color','r');
    hold on
    line2(:,1)=nset2(i,1);
    line2(:,2)=nset2(i,2);
    line2(:,3)=0;
    line2(:,4)=nset2(i,3);
    line2(:,5)=nset2(i,4);
    line2(:,6)=0; 
    write_dxf_line(fid1,line2);
end
write_dxf_end(fid1);


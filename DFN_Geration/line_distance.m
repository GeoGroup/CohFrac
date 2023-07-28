function mina = line_distance(points1,points2)
g=points1(2)-points1(4);
h=points1(1)-points1(3);
if g<0|h<0
    t1=-1:0.2:0;
else
    t1=0:0.2:1;
end
x1=points1(1)+h*t1;
y1=points1(2)+g*t1;
 
 
g1=points2(2)-points2(4);
h1=points2(1)-points2(3);
if g1<0 | h1<0
    t2=-1:0.1:0;
else
    t2=0:0.1:1;
end
x2=points2(1)+h1*t2;
y2=points2(2)+g1*t2;
 
n=length(t1);
m=length(t2);
 
juli=zeros(n,m);
for i=1:n
    juli(i,:)=sqrt((x2-x1(i)).^2+(y2-y1(i)).^2);%
end%获得每两个点之间的距离
[mina,mini]=min(juli(:));%找出距离中的最小值，及其单下标

end
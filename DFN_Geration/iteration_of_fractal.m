function [xpp,ypp,m]=iteration_of_fractal(x,y,d,L)
% 节理分形迭代
% L为迭代次数
N=length(x);
k=1;    %迭代第1次
for j=1:length(x)-1
    a(j)=(x(j+1)-x(j))/(x(N)-x(1));
    c(j)=((y(j+1)-y(j))-d*(y(N)-y(1)))/(x(N)-x(1));
    e(j)=(x(N)*x(j)-x(1)*x(j+1))/(x(N)-x(1));
    f(j)=(x(N)*y(j)-x(1)*y(j+1)-d*(x(N)*y(1)-x(1)*y(N)))/(x(N)-x(1));
    if j==1            
        for i=1:length(x)
            xp(k)=a(j)*x(i)+e(j);
            yp(k)=c(j)*x(i)+d*y(i)+f(j);
            k=k+1;
        end
    else
        for i=2:length(x)
            xp(k)=a(j)*x(i)+e(j);
            yp(k)=c(j)*x(i)+d*y(i)+f(j);
            k=k+1;
        end
    end
end

for k=1:L-1 % 迭代剩下的（L-1）次
    m=1;    
    for j=1:length(x)-1
        if j==1
            for i=1:length(xp)
                xpp(m)=a(j)*xp(i)+e(j);
                ypp(m)=c(j)*xp(i)+d*yp(i)+f(j);
                m=m+1;
            end
        else
            for i=2:length(xp)
                xpp(m)=a(j)*xp(i)+e(j);
                ypp(m)=c(j)*xp(i)+d*yp(i)+f(j);
                m=m+1;
            end
        end
    end
    xp=xpp;
    yp=ypp;
end
end
function visual_poly(poly)
%Visual polygons using a set points
poly(size(poly,1)+1,:)=poly(1,:);
plot(poly(:,1),poly(:,2),'-'); hold on;
end


n=100;
pts = rand(n,2);                                                        % locations~ U(0,1)
for i=1:n-1
    while norm(pts(i,:)-pts(i+1,:)) < 0.01
        pts(i,:) = rand(1,2);
    end
end
clear;
N=1e7;
A = [-.25,0];
B = [-A(1),A(2)];
C = [A(2),B(1)];

tri_points = [A(1),B(1),C(1);A(2),B(2),C(2)];

points = zeros(N,2);
points(1,1) = A(1) + (B(1)-A(1))*rand(1);
points(1,2) = C(1) + (C(2)-abs(points(1,1)))*rand(1);


for i = 1:N-1
    point = randi([1,3],1);
    points(i+1,:) = (points(i,:)' + tri_points(:,point)).'/2;
end

figure()
hold on
scatter(points(:,1),points(:,2), 5,'o','filled')
plot([tri_points(1,:),A(1)], [tri_points(2,:),A(2)])
%axis([min(tri_points(1,:)) max(tri_points(1,:)) min(tri_points(2,:)) max(tri_points(2,:))])
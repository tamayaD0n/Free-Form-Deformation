
[x, y, z] = sphere(10); % create a sphere with 10 by 10 surface

% make sphere inside the unit cubic of [0,1]*[0,1]*[0,1] rotate pi/4 clockwise around the z-axis
% sphere radius = 0.45, sphere center = (1/sqrt(2), 0, 0.5)
x = 0.45 * x + 0.5*sqrt(2);
y = 0.45 * y + 0;
z = 0.45 * z + 0.5;

% surface plot
figure;
surf(x, y, z);

axis equal; 
xlim([0 2]);
ylim([-1 1]);
zlim([0 1]);
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;


% turn n by n matrix x,y,z into vector x1,y1,z1 
x1 = zeros(121,1);
y1 = zeros(121,1);
z1 = zeros(121,1);

for a = 1:11
    for i = 1:11
        num = (a-1)*11+i;
        x1(num) = x(a,i);
        y1(num) = y(a,i);
        z1(num) = z(a,i);
    end
end

figure;
plot3(x1, y1, z1,'o');

%store x1,y1,z1 in sphere121.mat
save('sphere121.mat', 'x1', 'y1', 'z1');
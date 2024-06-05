%% https://dl.acm.org/doi/pdf/10.1145/15886.15903

%load points of surface
load('sphere121.mat', 'x1', 'y1', 'z1');
X_points = [x1,y1,z1];

%define origin point，unit vector of x,y,z axis
X_0 = [0,0,0];
X = [1,0,0];
Y = [0,1,0];
Z = [0,0,1];


%unit vector of s,t,u axis
S = [0.5*sqrt(2),-0.5*sqrt(2),0];
T = [0.5*sqrt(2),0.5*sqrt(2),0];
U = [0,0,1];

l = 1;
m = 2;
n = 3;

%%
P_ijk = zeros(l+1, m+1, n+1, 3);
for i = 0:l
    for j = 0:m
        for k = 0:n
            P_ijk(i+1,j+1,k+1,:) = X_0 + (i/l)*S + (j/m)*T + (k/n)*U;
        end
    end
end

P_ijk_o = P_ijk; %P_ijk_o grid of control points

P_ijk(2,1,1,2) = -0.5; %P_ijk grid of control points after deformation
P_ijk(2,3,1,2) = 1.5;

%% find X's coordinate in stu

X_mid = zeros(size(X_points,1),size(X_points,2));
for i = 1:size(X_points,1)
    [s,t,u] = convertToSTU(X_points(i,:),X_0,S,T,U);
    X_mid(i,:) = [s,t,u];
end


%% calculate X's coordiante in xyz after ffd
X_ffd = zeros(size(X_points,1),size(X_points,2));
for i = 1:size(X_points,1)
    Xffd = calculateXffd(l, m, n, X_mid(i,1), X_mid(i,2), X_mid(i,3), P_ijk);
    X_ffd(i,:) = Xffd;
end

%%
plot_x = zeros(11,11);
plot_y = zeros(11,11);
plot_z = zeros(11,11);

for a = 1:11
    for i = 1:11
        num = (a-1)*11+i;
        plot_x(a,i) = X_points(num,1);
        plot_y(a,i) = X_points(num,2);
        plot_z(a,i) = X_points(num,3);
    end
end
figure();
plotPijk(P_ijk_o);
surf(plot_x,plot_y,plot_z);
axis equal; 
xlim([-1 2]);
xlim([0 2]);
ylim([-1 1]);
zlim([0 1]);
grid on;


for a = 1:11
    for i = 1:11
        num = (a-1)*11+i;
        plot_x(a,i) = X_ffd(num,1);
        plot_y(a,i) = X_ffd(num,2);
        plot_z(a,i) = X_ffd(num,3);
    end
end
figure();
plotPijk(P_ijk);
surf(plot_x,plot_y,plot_z);
axis equal; 
xlim([0 2]);
ylim([-1 1]);
zlim([0 1]);
grid on;
    
    
%% function block
function [s,t,u] = convertToSTU(X,X0,S,T,U)
    TU = cross(T,U);
    numerator_S = dot(TU, (X - X0));
    denominator_s = dot(TU, S);

    s = numerator_S / denominator_s; 
    SU = cross(S,U);
    numerator_S = dot(SU, (X - X0));
    denominator_s = dot(SU, T);
    t = numerator_S / denominator_s; 
    
    ST = cross(S,T);
    numerator_S = dot(ST, (X - X0));
    denominator_s = dot(ST, U);
    u = numerator_S / denominator_s; 
    
    % or in short, 
    % note that in this case, X should be a 3*n matrix
    % or say point X expressed by column vector
    %     A = [S, T, U];
    %     P = X - X0;
    %     V = A \ P; %[s,t,u]=V(1),V(2),V(3);
    
end

function Xffd = calculateXffd(l, m, n, s, t, u, P)
    Xffd = zeros(1, 3); 
    for i = 0:l
        term1 = nchoosek(l, i) * (1 - s)^(l - i) * s^i;
        for j = 0:m
            term2 = nchoosek(m, j) * (1 - t)^(m - j) * t^j;
            for k = 0:n
                term3 = nchoosek(n, k) * (1 - u)^(n - k) * u^k;
                P_mid = [P(i+1,j+1,k+1,1),P(i+1,j+1,k+1,2),P(i+1,j+1,k+1,3)];
                Xffd = Xffd + term1 * term2 * term3 * P_mid;
            end
        end
    end
end

function plotPijk(P)
    for i = 1:size(P,1)
        for j = 1:size(P,2)
            for k = 1:size(P,3)
                plot3(P(i,j,k,1),P(i,j,k,2),P(i,j,k,3),'o','Color','b');
                hold on;
            end
        end
    end
end



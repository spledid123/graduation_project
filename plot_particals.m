%   根据粒子位置，画粒子
filename = 'data\bulk_pore_pm17_R_500_N_1000000_dT_1_2000\rxT_circle_T_250_.txt';
%   不同粒子的位置
A = readtable(filename);
len = size(A);
len = len(1);
%   粒子数
path_start = 1;
path_end = len;
rx = A.rx;
ry = A.ry;


%%   全画
figure;
axis equal;
hold on;
plot(rx,ry,'b.');
%%   可能存在不需要画的粒子
figure;
axis equal;
hold on;
t = (rx.^2 + ry.^2 - 500 ^ 2) < -1e-3;% 不用画的为0
x = rx .* t;
y = ry .* t;
% 绘制粒子图
plot(x,y,'b.');
plot(0,0,'w.');


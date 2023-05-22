%   主函数，描述一个粒子的轨迹
%   给定温度（决定速度），给定初位置x0,初速度gx0,V0,给定限定条件时间间隔N_T与sum_j
%   得到初始时间T=0到T=sum_j*N_T时间内，每次碰撞的位置表[时间 坐标x 坐标y 与边界的距离]([T r(1) r(2) sd])，
%   表的第一行是初始位置。
%   给出多孔介质的形状参数R,含义取决于scene.m
clc;
clear;

%   输入参数
N_T = 1e-1;         %   时间间隔
sum_j = 2e3;        %   粒子跑多少个时间间隔停止。
R = [1.1 1.1];      %   结构参数，取决于scene.m
Tem = 40;           %   绝对温度

x0 = [0 0];         %   给定初位置，在原点
gx0 = [];           %   初速度方向单位向量
V0 = 0;             %   初速度

x = [];
theta = 0;          %   速度方向角度
gx = [];
V = 0;              %   位置，速度方向

T = 0;              %   时间
T_ad = 0;           %   剩余吸附时间

%   输出单粒子轨迹表
filename = 'data\onepartical_long_path_cir\\rxT_circle_onepartical_sixcir2.txt';
varNames = {'T','rx','ry','sd'};
writetable(table(T,r(1),r(2),1,'VariableNames',varNames),filename,'WriteMode','append');

%   j决定粒子跑的时间，即：跑多少时间间隔
for j = 0:sum_j
    disp(j);        %   输出j,当进度条用。
    T1 = N_T;       %   onegas的参量
    if j == 0       %   T=0，初始化粒子的位置，速度，剩余吸附时间
        %   粒子的初始信息，给出初速度
        theta = rand() * 2 * pi;
        gx0 = [cos(theta) sin(theta)];
        V0 = Boltzmann(Tem);

        T_ad = 0;
        x = x0;
        gx = gx0;
        V = V0;
        continue;
    end

    [s,x,gx,V,T_ad] = onegas(x, x0, gx, V, T1, T_ad, Tem, filename, T, R);
    %    输入T时的位置与速度，给出T1后的位置与速度
    T = T + T1;
end





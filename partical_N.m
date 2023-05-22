%   主函数，描述N个粒子的运动，赋予速度，固定时间间隔记录粒子位置，速度
%   并行计算
%   考虑周期性边界
clc;
clear;

EPSILON = 1e-6;
EPSILON_SD = 1e-7;

R = [500 1.5 100];     %   几何参数
% 九宫格多孔介质参数R=[右边的圆半径 右边的圆半径+喉道宽度一半 右边的一排/列圆数量 左边的圆半径 左边的圆半径+喉道宽度一半 左边的一排/列圆数量]
%  R = [r,Rp,nump],r为bulk圆半径,圆心在原点,Rp为阵列的孔喉比,nump为圆阵列每行个数，圆的半径为1

N = 50000;            %   粒子个数
N_T = 1e-1;         %   时间间隔

Tem = 40;           %   绝对温度
theta = 0;          %   速度方向角度
x0 = [0 0];         %   原点，初位置
gx0 = [0 0];        %   初速度方向单位向量
V0 = 0;             %   初速度
T_ad = 0;           %   停留时间
T = 0;              %   记录时间

x = [];
gx = [];
V = 0;              %   位置，速度方向

d_x = ones(N,2);
d_gx = ones(N,2);
d_v = ones(N,1);
d_T = ones(N,1);    %   每个粒子的位置，速度，停留时间
d_sd = zeros(N,1);  %   sdf
dc_x = cell(N,1);   %   碰撞的位置
dc_y = cell(N,1);   
dc_t = zeros(N,1);  %   碰撞的次数

filenamesta = 'data\';
filenameend = '_.txt';
%  直接初始化
% for i = 1:N
%     %   角度均匀发射
%     theta = 2*pi*rand();
%     gx0 = [cos(theta) sin(theta)];
%     V0 = Boltzmann(Tem);
%     T_ad = 0;
%     d_x(i,:) = x0;
%     d_v(i) = V0;
%     d_gx(i,:) = gx0;
%     d_T(i) = T_ad;
% end
%   初始化，已经跑了一段时间，直接读取文件
% T = 300;    %   起始时刻
% filenamemid = num2str(T);
% filename = strcat(filenamesta,filenamemid,filenameend);
% A = readtable(filename);
% d_x = table2array(A(:,1:2));
% d_gx = table2array(A(:,3:4));
% d_v = table2array(A(:,5));
% d_T = table2array(A(:,6));
%   初始化，针对圆bulk内部有多孔介质的体系,向园内射入粒子，方向满足余弦，速度满足2D泄流
% r = R(1);
% for i = 1:N
%     %   角度均匀发射
%     theta = 2*pi*rand();
%     theta_v = theta + (asin(2*rand()-1));
%     gx0 = [cos(theta_v) sin(theta_v)];
%     V0 = Boltzmann(Tem);
%     T_ad = (i-1+250) * (1/200);   %  轮流发射
%     d_x(i,:) = [-r * cos(theta) -r * sin(theta)];
% %     d_x(i,:) = x0;
%     d_v(i) = V0;
%     d_gx(i,:) = gx0;
%     d_T(i) = T_ad;
%     d_sd(i) = scene(d_x(i,:),R);
% end

sumj = 5500 - 3000;
for j = 1:sumj %    时间循环
    parfor i = 1:N
        x = d_x(i,:);
        gx = d_gx(i,:);
        V = d_v(i);
        T_ad = d_T(i);
        sd = d_sd(i);
        if(~isnan(sd))
            [x, gx, V, T_ad, xcoll, tcoll, sd] = Ngas(R, x, gx, V, N_T, T_ad, Tem);
        end
%         %   考虑周期性边界，归一
%         t = 2 * R(3) * R(2);           %   普通区块
%         t = 2 * R(3) * R(2) + 2 * R(7);%    考虑有空余边界的区块
%         t = 40 * 1.5;
% %         tx = 2 * numx * R(1);
% %         ty = 2 * numy * R(2);
%         tx = t;
%         ty = t;
%         kx = floor(x(1)/tx);
%         ky = floor(x(2)/ty);
%         if(mod(kx,2) == mod(ky,2))    %   x > 0
%             x(1) = x(1) - kx * tx;
%             x(2) = x(2) - ky * ty;
%         else                          %   x < 0
%             x(1) = x(1) - (kx + 1) * tx;
%             x(2) = x(2) - ky * ty;
%         end
        %   对碰撞点的处理,记录下来
%         for m = 1:tcoll
%         kx = floor(xcoll(m,1)/tx);
%         ky = floor(xcoll(m,2)/ty);
%         %   归一
%         if(mod(kx,2) == mod(ky,2))
%             xcoll(m,1) = xcoll(m,1) - kx * tx;
%             xcoll(m,2) = xcoll(m,2) - ky * ty;
%         else  
%             xcoll(m,1) = xcoll(m,1) - (kx + 1) * tx;
%             xcoll(m,2) = xcoll(m,2) - ky * ty;
%         end
%         end
%         %   临时功能：碰撞点输入文件
%         if(tcoll > 0)
%         dc_x{i} = [dc_x{i};xcoll(:,1)];
%         dc_y{i} = [dc_y{i};xcoll(:,2)];
%         dc_t(i) = dc_t(i) + tcoll;
%         end
        %
        d_x(i,:) = x;
        d_v(i) = V;
        d_gx(i,:) = gx;
        d_T(i) = T_ad;
        d_sd(i) = sd;
    end
    %  记录粒子位置，速度
    if(mod(j,10)==0)
        fprintf("进度%f%%\n",j*100/sumj);
    end
    T = T + N_T;
    filenamemid = num2str(T);
    filename = strcat(filenamesta,filenamemid,filenameend);
    varNames = {'rx','ry','gx','gy','v','Tad','sd'};
    writetable(table(d_x(:,1),d_x(:,2),d_gx(:,1),d_gx(:,2),d_v,d_T,d_sd,'VariableNames',varNames),filename,'WriteMode','append');
end

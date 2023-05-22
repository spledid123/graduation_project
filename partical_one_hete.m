%   圆半径非均匀的阵列，粒子在其中运动，形成轨迹表
%   半径实时生成，半径均匀分布
%   给定最大碰撞次数
clc;
clear;
x0 = [0 0];%    初始位置
gx0 = [0 0];%   初始速度方向
N_pmax = 2e4;%  最大碰撞次数
global Get_ran_r_flag;% 是否生成过半径
global Get_ran_r_r;%    存储所有半径数据
global Get_ran_r_maxflag;%  Get_ran_r_r的大小
Get_ran_r_flag = [];
Get_ran_r_r = [];
Get_ran_r_maxflag = 0;

R = [1.2 1.2];%   x,y方向圆心距离
Rmax = [1 1.1 1.05 1.15];%    最大半径 
Rmin = [1 0.9 0.95 0.85];%  最小半径

s = [];%    经过一定次数的碰撞后，粒子的位移大小
x = [];%    经过一定次数的碰撞后，粒子的x位移
y = [];%    经过一定次数的碰撞后，粒子的y位移
gx = [];

num_par = 4;
%   文件名
filenameend = '_.txt';
filenamesta = 'data\2\rxT_circle_onepartical_sixcir_Rmin_';
filemid1 = "_Rmin_";

for m = 1:num_par
    filenamemid1 = num2str(Rmax(m));
    filenamemid2 = num2str(Rmin(m));
    file3 = num2str(m);
    filename = strcat(filenamesta,filenamemid1,filemid1,filenamemid2,file3,filenameend);
    theta = pi/4;
    gx0 = [cos(theta) sin(theta)];
    [filename, rx, gx, s(m)] = fun_partical_one_hete(x0, gx0, N_pmax, filename, R, Rmax(m), Rmin(m));
    x(m) = rx(1);
    y(m) = rx(2);
end

function [filename, x, gx, s] = fun_partical_one_hete(x0, gx0, N_pmax, filename, R, Rmax, Rmin)
%   这是一个函数，描述一个粒子在非均匀圆阵列的轨迹，方便并行
%   我只需要这个粒子每次碰撞的位置
%   输入参数：[x0,gx0,N_pmax,filename,R],分别为粒子的初始位置，运动的初始方向,最大碰撞次数(判停条件)，输出文件名,形状参数
%   输出参数:[filename,x,gx,s],表示文件名,末位置,末运动方向,末运动位置与原点的距离
%   Rmax,Rmin为最大最小半径

EPSILON = 1e-6;%    小量
EPSILON_SD = 1e-7;
LOOP_MAX = 1e7;%    防止无限循环

x = x0;
x1 = x;
gx = gx0;
s = 0;
sd = 0;

T = 0;% 没用的变量，占位；
N_p = 0;%   碰撞次数

varNames = {'T','rx','ry','sd'};
writetable(table(T,x(1),x(2),sd,'VariableNames',varNames),filename,'WriteMode','append');%    第一行的位置为初位置

i = 0;
while(N_p < N_pmax && i < LOOP_MAX)
    i = i + 1;
    sd = scene_hete(x, R, Rmax, Rmin);
    if(sd + EPSILON_SD <= 0)% 卡进固体里去了
        error("sd < 0\n,现位置x=%f,y=%f;上次运动开始时的位置x=%f,y=%f,sd=%f",x(1),x(2),x1(1),x1(2),sd);
    elseif(sd-EPSILON < 0)          % 在边界
        writetable(table(T,x(1),x(2),scene_hete(x,R,Rmax,Rmin)),filename,'WriteMode','append');
        N_p = N_p + 1;
        n = gradient_hete(x, R,Rmax,Rmin);
        gx = reflect(gx, n);
        x = x + gx .* 2 * EPSILON;% 避免陷在边界
        s = x - x0;
        s = (s(1)^2+s(2)^2)^0.5;
        continue;
    else
        %   非特殊情况，不在边界附近
        x1 = x;
        x = x + gx * sd;
        s = x - x0;
        s = (s(1)^2+s(2)^2)^0.5;
    end
end
if(i == LOOP_MAX)
    sprintf('LOOP_MAX\n');
end
end



function sd = scene_hete(x,R,Rmax,Rmin)
%   返回点到边界的最短距离
%   x=[x,y]为坐标
%   R=[Rx,Ry]为多孔介质形状参数
%   Rmax,Rmin为get_ran_r函数的参数

% %   非均匀圆形阵列
%   圆的半径为1,于是x方向孔的半径为Rx-1;y方向孔的半径为Ry-1;
%   于是x方向圆心坐标为(Rx)*(2*kx+1);y方向圆心坐标为Ry*(2*ky+1);
kx = floor((x(1)/R(1)-1)/2);
ky = floor((x(2)/R(2)-1)/2);

rxmin = R(1) * (2 * kx + 1);
rxmax = R(1) * (2 * kx + 3);
rymin = R(2) * (2 * ky + 1);
rymax = R(2) * (2 * ky + 3);

a=[];
a(1) = circleSDF(x,[rxmin,rymin],get_ran_r(kx,ky,R,Rmax,Rmin));
a(2) = circleSDF(x,[rxmax,rymin],get_ran_r(kx+1,ky,R,Rmax,Rmin));
a(3) = circleSDF(x,[rxmin,rymax],get_ran_r(kx,ky+1,R,Rmax,Rmin));
a(4) = circleSDF(x,[rxmax,rymax],get_ran_r(kx+1,ky+1,R,Rmax,Rmin));

sd = min(a);

end
function n = gradient_hete(x, R, Rmax, Rmin)
%   求边界的法向量n=(nx,ny)
%   x=[x,y]为边界点的坐标，R为形状参数，为函数scene的参数
EPSILON = 1e-6;

n(1) = (scene_hete([x(1) + EPSILON x(2)],R, Rmax,Rmin) - scene_hete([x(1) - EPSILON x(2)],R,Rmax,Rmin)) * 0.5 / EPSILON;
n(2) = (scene_hete([x(1) x(2) + EPSILON],R,Rmax,Rmin) - scene_hete([x(1) x(2) - EPSILON],R,Rmax,Rmin)) * 0.5 / EPSILON;

len = sqrt(n(1)^2+n(2)^2);
n = n ./ len;
end
function sd  = get_ran_r(kx, ky, R, Rmax, Rmin)
%   生成随机分布的圆阵列，即生成圆半径
%   R表示圆心距离
%   半径均匀分布
%   Rmax,Rmin表示半径最大最小值。
%   已经形成的记录下来,分别用全局变量Get_ran_r_flag表示(kx,ky)是否生成过，1为生成过
%   Get_ran_r_maxflag为Get_ran_r_flag的大小
%   全局变量Get_ran_r_r(kx,ky)表示生成的变量大小
global Get_ran_r_flag;
global Get_ran_r_r;
global Get_ran_r_maxflag;
% Get_ran_r_flag=[],记录方式为：(0,0),(1,0),(0,1),(-1,0),(0,-1),(2,0),(1,1),(0,2),(-1,-1),(-2,0),...
t = fun_fromktoflag(kx,ky);
if(t > Get_ran_r_maxflag)
    for i = Get_ran_r_maxflag + 1:t
        Get_ran_r_flag(i) = 0;
    end
    Get_ran_r_maxflag = t;
end
if(Get_ran_r_flag(t) == 0)
    sd = rand()*(Rmax - Rmin) + Rmin;
    Get_ran_r_r(t) = sd;
    Get_ran_r_flag(t) = 1;
    
else
    sd = Get_ran_r_r(t);
end
end
function x = fun_fromktoflag(kx,ky)
    %得到Get_ran_r_flag的序号
    if(kx==0 && ky==0)
        x = 1;
        return;
    end
    t = abs(kx)+abs(ky);
    x = 2*t*(t-1)+1;
    if(kx == 0)
        if(ky > 0)
            x = x + ky + 1;
        else
            x = x - ky * 3 + 1;
        end
    elseif(ky == 0)
        if(kx > 0)
            x = x + 1;
        else
            x = x - kx * 2 + 1;
        end
    else
        if(kx > 0 && ky > 0)
            x = x + ky + 1;
        elseif(kx < 0 && ky > 0)
            x = x + t - kx + 1;
        elseif(kx < 0 && ky < 0)
            x = x + t * 2 - ky + 1;
        else
            x = x + t * 3 + kx + 1;
        end
    end
end
function [kx,ky] = fun_fromflagtok(x)
if(x == 1)
    kx = 0;
    ky = 0;
    return;
end
t = (-2+(4-8*(1-x))^0.5)/4;
t = floor(t);
if(t < 0)
    error(0);
end
a = 2*t*(t-1)+1;
a = x-a;
if(a == 0)
    ky = 0;
    kx = t;
    return;
end
t = t + 1;
kt = floor(a/t);
a = a - kt * t;
if(kt == 0)
    ky = a - 1;
    kx = t - kx;
elseif(kt == 1)
    kx = -(a - 1);
    ky = t + kx;
elseif(kt == 2)
    ky = -(a-1);
    kx = -(t+ky);
else
    kx = a-1;
    ky = -(t-kx);
end
end
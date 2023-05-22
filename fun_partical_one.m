function [filename, x, gx, s] = fun_partical_one(x0, gx0, N_pmax, filename, R)
%   这是一个函数，描述一个粒子的轨迹，方便并行
%   我只需要这个粒子每次碰撞的位置
%   输入参数：[x0,gx0,N_pmax,filename,R,Tem],分别为粒子的初始位置，运动的初始方向,最大碰撞次数(判停条件)，输出文件名,形状参数,温度(速度参数)
%   输出参数:[filename,x,gx,s],表示文件名,末位置,末运动方向,末运动位置与原点的距离

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
    sd = scene(x, R);
    if(sd + EPSILON_SD <= 0)% 卡进固体里去了
        error("sd < 0\n,现位置x=%f,y=%f;上次运动开始时的位置x=%f,y=%f,sd=%f",x(1),x(2),x1(1),x1(2),sd);
    elseif(sd-EPSILON < 0)          % 在边界
        writetable(table(T,x(1),x(2),scene(x,R)),filename,'WriteMode','Append', 'WriteVariableNames',false,'WriteRowNames',true);
        N_p = N_p + 1;
        n = gradient(x, R);
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
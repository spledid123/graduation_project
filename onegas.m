function [s, x, gx, V, T_newad] = onegas(x0, x00, gx0, V0, T0, T_ad, Tem, filename, Tall, R)
%   给定时间间隔T0,给定粒子的剩余吸附时间T_ad,给定时间为Tall的位置x0,速度gx0,V0,初始位置x00
%   给定绝对温度Tem
%   给定结构参数R
%   输出时间为Tall+T时，粒子的末位置x,gx,与初始位置x00的距离s，速度gx,V,以及剩余吸附时间T_newad
%   每次碰撞，将对应的时间，位置计入文件filename中,形成轨迹表
EPSILON = 1e-6;     %   小量
EPSILON_SD = 1e-7;
LOOP_MAX = 1e7;     %   最大循环次数

sd = 0;             %   与边界的距离
n = [0 0];          %   边界的法线方向
T = 0;              %   剩余需要考虑的时间。
T1 = 0;
Ta = Tall;

s = sqrt((x0(1) - x00(1))^2+(x0(2) - x00(2))^2);
x = x0;
gx = gx0;
V = V0;
x1 = [];              %   上一次移动前

%   若剩余吸附时间大于时间间隔T，位置不变
if T_ad > T0
    T_newad = T_ad -T;
    return
else    %   脱附
    T = T0 - T_ad;
    T_newad = 0;
end

i = 0;
while(T > EPSILON && i < LOOP_MAX)
    i = i + 1;
    sd = scene(x, R);           %   与边界的距离
    if(sd + EPSILON_SD <= 0)    %   越过边界
        sprintf('endT = %f\n',T);
        error("现位置x=%f,y=%f;上次运动开始时的位置x=%f,y=%f,sd=%f",x(1),x(2),x1(1),x1(2),sd);
    elseif(sd - EPSILON < 0)          %     在边界
        Ta = Tall + T0 - T;       %     在边界脱附时刻的时间
        writetable(table(Ta,x(1),x(2),scene(x,R)),filename,'WriteMode','Append');
        n = gradient(x, R);
        gx = reflect(gx, n);      %     脱附的速度
        V = Boltzmann(Tem);
        x = x + gx .* 2 * EPSILON; %    先移动一段距离，避免在边界移动太慢或者跨过边界
        %   乘2是为了保证x不在边界，于是再次进入循环后不会产生二次吸附
        T = T - 2 * EPSILON / V;
        T_newad = T_absorb();       %   吸附，给吸附时间
        if T_newad > T
            T_newad = T_newad - T;
            s = sqrt((x(1) - x00(1))^2+(x(2) - x00(2))^2);
            return
        else
            T = T - T_newad;
            T_newad = 0;
        end
        s = sqrt((x(1) - x00(1))^2+(x(2) - x00(2))^2);
        continue;
    else    %   在边界外
        x1 = x;
        T1 = T - sd / V;
        if(T1 * V > EPSILON)    %   往前sd后依然可以继续
            x = x + gx * sd;
            T = T - sd / V;
        else                    %   剩余时间不够，不能向前sd的距离
            sd = V * T;
            x = x + gx * sd;
            T = 0;
        end
    end
end
if(i == LOOP_MAX)   %   死循环
    sprintf('LOOP_MAX\n');
end

s = sqrt((x(1) - x00(1))^2+(x(2) - x00(2))^2);
end

























function [x, gx, V, T_newad, xcoll,tcoll, sd] = Ngas(R, x0, gx0, V0, T, T_ad, Tem, varargin)
%   给定一定的时间间隔，粒子的原位置，速度,输出粒子的末位置，速度
%   考虑吸附
%   输入：R为几何参数,x0为原位置，gx0，V0为原速度,T为时间间隔,T_ad为剩余吸附时间,Tem为温度
%   输出：x为末位置,gx，V为末速度,T_newad为末剩余吸附时间
%   xcoll为所有碰撞的点,tcoll为碰撞次数
%   pm为可选参数，多孔介质信息
p = inputParser;            % 函数的输入解析器
addOptional(p,'pm',0);
parse(p,varargin{:});

EPSILON = 1e-6;%    时间小量
EPSILON_SD = 1e-1;% 距离小量,非像素计算可以是1e-6
LOOP_MAX = 1e7;

sd = 0;%    到边界的距离
n = [0 0];% 边界法线

x = x0;
gx = gx0;
V = V0;
x1 = x0;   %   记录上一次移动的起点
xcoll = [];
tcoll = 0;


if T_ad > T
    T_newad = T_ad -T;
    return
else
    T = T - T_ad;
    T_newad = 0;
end

i = 0;
while(T > EPSILON && i < LOOP_MAX)
    i = i + 1;
    sd = scene(x, R, p.Results.pm);
    if(sd + EPSILON_SD <= 0)
        error("out of bound\n位置x=%f,y=%f,sd=%f\n速度方向gx(1)=%f,gx(2)=%f\n上次移动的位置x1=%f,y1=%f,sd=%f",x(1),x(2),sd,gx(1),gx(2),x1(1),x1(2),scene(x1,R,p.Results.pm));
    elseif(sd - EPSILON_SD < 0)          % 在边界
        tcoll = tcoll + 1;% 记录碰撞点
        xcoll(tcoll,:) = x;
        n = gradient(x, R, EPSILON_SD * 10, p.Results.pm);
        gx = reflect(gx, n);
        V = Boltzmann(Tem);
        x = x + gx * 2 * EPSILON_SD;
        T = T - 2 * EPSILON_SD / V;
        T_newad = T_absorb();
        if T_newad > T
            T_newad = T_newad -T;
            return;
        else
            T = T - T_newad;
            T_newad = 0;
        end
        continue;
    else
        x1 = x;
        T1 = T - sd / V;
        if(T * V - sd > -EPSILON_SD)
            sd = sd - EPSILON_SD;
            x = x + gx * sd;
            T = T - sd / V;
        else
            sd = V * T;
            x = x + gx * sd;
            T = 0;
        end
    end
end
if(i == LOOP_MAX)
    sprintf('LOOP_MAX\n');
end
sd = scene(x, R, p.Results.pm);
if(sd + EPSILON_SD <= 0)
    error("out of bound\n位置x=%f,y=%f,sd=%f\n速度方向gx(1)=%f,gx(2)=%f\n上次移动的位置x1=%f,y1=%f,sd=%f",x(1),x(2),sd,gx(1),gx(2),x1(1),x1(2),scene(x1,R,p.Results.pm));
end
end

function [isout, time] = onegaspath(x, V, gx, ppm, fsn)
%   给定某个粒子的位置x,速率V,速度方向单位向量gx,多孔介质sdf信息ppm,判定是否在多孔介质内部函数fsn
%   isout表示粒子是否直接离开，1表示不需要碰撞；time表示花费时间
%   fsn(x,y)求边界的sdf,在里面是负数，在外面是正数
%   若粒子在外面，则输出isout = 0;
%   可以通过改变fsn调换"内外"
time = 0;
EPSILON_SD = 1e-1;
LOOP_MAX = 1e7;%    防止无限循环

if(fsn(x(1), x(2)) > 0)
    isout = 0;
    time = 0;
else
    i = 0;
    while(i < LOOP_MAX)
        i = i + 1;
        sd = scene(x, R, ppm);
        if(sd + EPSILON_SD <= 0)% 卡进固体里去了
            error("sd < 0\n,现位置x=%f,y=%f;上次运动开始时的位置x=%f,y=%f,sd=%f",x(1),x(2),x1(1),x1(2),sd);
        elseif(sd - EPSILON_SD < 0)          % 在边界
            isout = 0;
            break;
        else    %   非特殊情况，不在边界附近
            sd = sd - EPSILON_SD;
            x1 = x;
            x = x + gx * sd;
            if(fsn(x(1), x(2)) > 0)
                isout = 1;
                time = time + (-fsn(x1(1), x1(2))/(fsn(x(1), x(2)) - fsn(x1(1), x1(2)))) * sd / V;
            end
            time = time + sd/V;
        end
    end
if(i == LOOP_MAX)
    error('LOOP_MAX\n');
end
end
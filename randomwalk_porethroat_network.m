%   随机游走网格，每个格点是孔，而连接两个孔的线算作喉
%   每一步（单位时间），有一定概率从当前位置前往下一个临近位置，即若当前在孔，则可能留在原地，或者前往临近喉；
%   若当前在喉，则可能留在原地，或者前往临近孔；
%   当前在孔，某一步在原地的概率为Pp;当前在喉，某一步留在原地的概率为Pt;
%   在孔内的单步时间为Tp,在喉内的单步时间为Tt
%   关于计时：每一步以终点为准，若起点为孔，记为该孔时间加1；若起点为喉，临近两孔分别有一半概率时间加1。
%   每次改变计时的孔，时间另计：目的：统计单孔停留时间。
%   初始点在某个孔
%   通过随机游走模拟得到单孔停留时间的概率分布，单步时间可设为相同的或不同的
%   计算概率分布，得到pdf
%   蒙特卡洛模拟得到扩散系数

clc;
clear;
T = [];                     %   单孔停留时间
T0 = 0;                     %   本步时间
flag = 0;                   %   上一步在孔(0),在喉(1)
MAX_STEP = 1e7;             %   走的总步数
%   输入参数
Tp = 1.4334e-7;             %   在孔内的单步时间
Tt = 2.0702e-8;              %   在喉内的单步时间
pp = 0.3;           %   单次碰撞孔离开孔的概率
pt = 0.323;               %   单次碰撞喉离开喉的概率

%%   随机游走模拟，孔喉时间不同
T(1) = 0;
j = 1;                      %   计时的孔总数
for i = 1:MAX_STEP
    if(flag == 0)           %   上一步在孔
        T0 = Tp;            %   孔->孔
        r = rand();
        if(r < pp)      %   孔->临近喉
            flag = 1;
            T0 = Tt;
        end
        T(j) = T(j) + T0;
    else                    %   上一步在喉
        rs = rand();
        if(rs < 0.5)        %   一半概率进入下一个孔
            j = j + 1;
            T(j) = 0;
        end        
        T0 = Tt;            %   喉->喉
        r = rand();
        if(r < pt)          %   喉->临近孔
            flag = 0;
            T0 = Tp;
        end
        T(j) = T(j) + T0;
    end
end
%%  随机游走模拟，统计单元碰撞次数
N = [];
j = 1;                      %   计时的孔总数
N(1) = 0;
flag = 0;
for i = 1:MAX_STEP
    if(flag == 0)           %   上一步在孔
        r = rand();
        if(r < pp)      %   孔->临近喉
            flag = 1;
            N(j) = N(j) + 1;
        end
        N(j) = N(j) + 1;
    else                    %   上一步在喉
        rs = rand();
        if(rs < 0.5)        %   一半概率进入下一个孔
            j = j + 1;
            N(j)=1;
        end        
        N(j) = N(j)+1;            %   喉->喉
        r = rand();
        if(r < pt)          %   喉->临近孔
            flag = 0;
        end
    end
end
%%  随机游走模拟，计算扩散系数,没写
N = 1e4;    %   粒子数
T(1) = 0;
j = 1;                      %   计时的孔总数
for i = 1:MAX_STEP
    if(flag == 0)           %   上一步在孔
        T0 = Tp;            %   孔->孔
        r = rand();
        if(r < pp)      %   孔->临近喉
            flag = 1;
            T0 = Tt;
        end
        T(j) = T(j) + T0;
    else                    %   上一步在喉
        rs = rand();
        if(rs < 0.5)        %   一半概率进入下一个孔
            j = j + 1;
            T(j) = 0;
        end        
        T0 = Tt;            %   喉->喉
        r = rand();
        if(r < pt)          %   喉->临近孔
            flag = 0;
            T0 = Tp;
        end
        T(j) = T(j) + T0;
    end
end
%%   随机游走模拟，孔喉的单步时间相同，修正概率
T0 = 1.056e-7;               %   单步时间 
%   根据指数分布的特点修正概率
fp = -log(1-pp)/Tp;
ft = -log(1-pt)/Tt;
Pp = 1-exp(-fp*T0);
Pt = 1-exp(-ft*T0);

T(1) = 0;
for i = 1:MAX_STEP
    if(flag == 0)           %   上一步在孔
        r = rand();
        if(r > 1 - Pp)      %   孔->临近喉
            flag = 1;
        end
        T(j) = T(j) + T0;
    else                    %   上一步在喉
        rs = rand();
        if(rs < 0.5)        %   一半概率进入下一个孔
            j = j + 1;
            T(j) = 0;
        end        
        r = rand();
        if(r > Pt)          %   喉->临近孔
            flag = 0;
        end
        T(j) = T(j) + T0;
    end
end

%%   根据随机游走的结果画概率分布图
[xx,yy] = plot_distribution(T,0,50);
hold on;
figure;
plot(xx,yy);
hold on;
box on;
title("单孔停留时间分布","FontSize",14);
xlabel("单孔停留时间","FontSize",14);
ylabel("概率密度","FontSize",14);
%%   计算概率密度,即pdf，单步时间孔喉相等
Pp = pp;
Pt = pt;
Pt = 1 - Pt;    % 单步留在喉内的概率  
k = 1:70;%  在某个格点内的步数
tx = k * T0;
y = [];
for i=1:length(k)
    y(i) = P(k(i), Pp, Pt, 1, 1);
end
figure;
hold on;
plot(k,y);% 概率
set(gca,'yscale','log');
box on;
title("单孔碰撞次数分布","FontSize",14);
xlabel("单孔碰撞次数","FontSize",14);
ylabel("概率","FontSize",14);
y = y  / T0;%   概率密度
figure;
hold on;
plot(tx,y);
set(gca,'yscale','log');
box on;
title("单孔停留时间分布","FontSize",14);
xlabel("单孔停留时间","FontSize",14);
ylabel("概率密度","FontSize",14);
%%  展示概率分布的特点：双线性
Pt = 1-Pt;
for i=1:70
y2(i) = Pl(k(i), Pp, Pt);
end
y2 = y2 / T0;
y1 = y - y2;
figure;
hold on;
plot(x,y2,'r');
plot(x,y1,'b');
plot(x,y,'g');
box on;




function p = P4(k2, a, Pp, Pt)
p = 0;
c = 0;
for b = 0:k2 - 2 * a + 1
    c = k2 - 2 * a + 1 - b;
    p=p + nchoosek(b + a - 1, a - 1) * nchoosek(c + a - 1, a - 1) * (1 - Pp) ^ b * (Pt / 2)^c;
end
end

function p = P2(k2, Pp, Pt)
p=0;
for a = 1:floor((k2 + 1) / 2)
    p = p + Pp ^ a * ((1 - Pt)/2) ^ (a - 1) * P4(k2 , a, Pp, Pt);
end
end

function p = P1(k1, Pp, Pt)
p = (Pt / 2) ^ (k1 - 1) * (1 - Pt) / 2;
end

function p = Pr(k, Pp, Pt)
p = 0;
for k1 = 1:k - 1
    p = p + P1(k1, Pp, Pt) * P2(k - k1, Pp, Pt) * 1 / 2;
end
end

function p = Pl(k, Pp, Pt)
p = (Pt / 2) ^ k / 2;
end

function p=P(k, Pp, Pt, a, b)
p=a*Pl(k,Pp,Pt)+b*Pr(k,Pp,Pt);
end



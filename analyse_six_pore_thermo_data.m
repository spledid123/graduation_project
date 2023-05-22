%   分析在六边形孔喉结构多孔介质内气体与空腔内气体在热力学上的的区别
%   热力学性质由单粒子轨迹得到
%   以喉道中心处垂直线为分界线，分析单孔的平均数据,进而得到热力学数据
%   分析不同赋值不同路径给算法带来的误差
%   六边形孔喉结构的孔隙度

clc;
clear;
R = 1.1;                %   （喉道宽度+孔的直径）/孔的直径
a = 0.25;               %   喉道长度/孔的直径
sumn = 3;               %   路径文件数量
T_ab = 0;               %   平均吸附时间  
NUM_Partical = 100;     %   对每一条路径的赋值(速度)次数

meanDT_path = zeros(1,sumn);    %  单条路径的单孔平均时间
meanDV_path = zeros(1,sumn);    %  单条路径的平均越孔法向速率
meanDT = 0;                     %    单孔平均时间
meanDV = 0;                     %    平均越孔法向速率

% 文件名，从对应(R,a)取多条路径
filenamesta = "data\onepartical_long_path_six\rxT_circle_onepartical_sixcir_R_";
filenamemid = "_a_";
filenamemid2 = "(";
filenameend = ").txt";
filename1 = num2str(R);
filename2 = num2str(a);
%   为画图做准备，获取坐标区
figure;
hold on;
ax1 = gca;
figure;
hold on;
ax2 = gca;
figure;
hold on;
ax3 = gca;
figure;
hold on;
ax4 = gca;
for n = 1:sumn  %   不同的路径
    mDT = zeros(1,Npore);           %  单次赋值的前i个孔的单孔平均时间
    mDV = zeros(1,Npore);           %  单次赋值的前i个孔的平均越孔法向速率
    meanDT_one = zeros(1,Npore);    %  单条路径的前i个孔的单孔平均时间
    meanDV_one = zeros(1,Npore);    %  单条路径的前i个孔的平均越孔法向速率
    col = [rand() rand() rand()];   %  随机生成颜色,每一条轨迹一个颜色
    filename = strcat(filenamesta,filename1,filenamemid,filename2,filenamemid2,int2str(n),filenameend);
    %    返回每一次碰撞的参数
    [DR, kx, ky, TH, eX, eY] = per_walk(filename, [R, 1.1], a); %   修正方向角
    for i = 1:length(DR)                                        
        TH(i) = get_th(TH(i), kx(i), ky(i), kx(i+1), ky(i+1));
    end
    len = length(DR);
    %   考虑的路径范围
    start = 1;      
    send = len;
    for ii = 1:NUM_Partical %   同一条路径不同的速度
        %   得到孔为单位的数据
        [DT,V,Vth,Npore,meanlenDT] = times_one_particalpath(DR(start:send), kx(start:send+1), ky(start:send+1), TH(start:send), T_ab);
        V = V.*abs(cos(Vth));               %   越孔法向速度
        %   前i个孔的单孔平均时间
        for i=1:length(V)
            mDT(i) = sum(DT(1:i))/i;
        end
        %   前i个孔的平均越孔法向速率
        for i=1:length(V)
            mDV(i) = sum(V(1:i))/i;
        end
        %   对所有赋值求和
        meanDT_one = meanDT_one + mDT;
        meanDV_one = meanDV_one + mDV;
        %   画图，每一条线代表一次赋值
        plot(ax1,mDT,'Color',col);
        set(gca,'xscale','log');
        title('不同赋值下前i个孔的平均停留时间','FontSize',14);
        plot(ax2,mDV,'Color',col);
        set(gca,'xscale','log');
        title('不同赋值下前i个孔的平均越孔法向速度','FontSize',14);
    end
    %   一条路径在多次赋值的平均值
    meanDT_one = meanDT_one./NUM_Partical;
    meanDV_one = meanDV_one./NUM_Partical;
    %   画图，每一条线代表一条路径
    plot(ax3,meanDT_one,'Color',col);
    set(gca,'xscale','log');
    title('不同路径下前i个孔的平均停留时间','FontSize',14);
    plot(ax4,meanDV_one,'Color',col);
    set(gca,'xscale','log');
    title('不同路径下前i个孔的平均越孔法向速度','FontSize',14);
    %   该条路径所有孔的平均停留时间和法向速度
    meanDT_path(n) = meanDT_one(end);
    meanDV_path(n) = meanDV_one(end);
end
%   平均值
meanDT = sum(meanDT_path)/sumn;
meanDV = sum(meanDV_path)/sumn;
%   不同路径的误差
ef_DT = (max(meanDT_path)-min(meanDT_path))/meanDT;
ef_DV = (max(meanDV_path)-min(meanDV_path))/meanDT;
%%   画误差线表征不同路径的误差
plot(ax1,[1 Npore],[meanDT*0.9 meanDT*0.9],'b');
plot(ax3,[1 Npore],[meanDT*0.9 meanDT*0.9],'b');
plot(ax2,[1 Npore],[meanDV*1.1 meanDV*1.1],'b');
plot(ax4,[1 Npore],[meanDV*1.1 meanDV*1.1],'b');
%%  越孔速度方向与法相速率的分布
[xx,yy] = plot_distribution(Vth,0,100);
figure;
hold on;
plot(xx,yy,'r--','DisplayName','越孔速度方向分布');
title('越孔速度方向分布','fontsize',14);
set(gca,'xlim',[-pi/2 pi/2]);
xx = -pi/2:0.01:pi/2;
yy = cos(xx)/2;
plot(xx,yy,'b-','DisplayName','余弦分布');
%%      计算孔内分子数密度(压强)/空腔的值
r = 1;
t = (1 + a) * r;%   六边形边长的一半
h = 2*r*(R(1)-1);
len = 2*a*r;
A = pi + 1.5 * h * len;%    单孔面积
A = A * 1e-10;%    1m:10微米
h = h * 1e-5;

k = 1.380649e-23;
NA = 6.02e23;
m = 18e-3/NA;


%   三维玻尔兹曼分布(柱形)的比较
%   v = (8*k*40/pi/m)^0.5;
%   相同壁面碰撞频率下，孔内的实际分子数密度/对应的bulk分子数密度
%   设宽度为w,hw(n(bulk分子数密度)*v(平均速度)/4)=N(孔内分子数)/3(有3个开口)/meanDT;于是孔内分子数密度为N/(Aw)
%   f = 3*v*meanDT*h/4/A;

%   二维玻尔兹曼分布的比较
v = (pi*k*40/2/m)^0.5;
%   h(nv/pi)=N/3/meanDT,于是面积分子数密度为N/A,孔内分子数密度相对于外界bulk的数密度
f_2D = 3*v*meanDT*h/A/pi;
%    若是有吸附，需要修正，乘上（1-单步平均吸附时间-单步平均时间）
f_2D_ad = f_2D * (1 - T_ab/meanlenDT);
%  压强的比
fp_2D = 2 * meanDV * v / pi / (k * 40 / 2 / m);
%%  六边形孔喉结构的孔隙度
R = 1.1;
a = 0.25;
phi = get_poro_six(R, a);



function [DT, V, Vth, Npore, meanlenDT] = times_one_particalpath(DR, kx, ky, TH, T_ab)
%   输入某个粒子的多步路径长度DR,以及每次碰撞点的位置编号kx,ky(判断是否在同一个孔),
%   每次碰撞的修正后方向角TH，平均吸附时间T_ab
%   输出单孔时间DT与经过的孔的总数Npore
%   V为每个孔出去的速度大小。Vth为与法向的夹角
%   平均每一步meanlenDT的时间

%   赋予时间，没有吸附
filename_v = 'data\v_10000\Blotzmann_water_T_40_2D.txt';
Av = readtable(filename_v);
Av = table2array(Av);
len = length(DR) + 1;
T = zeros(1,len);           %   每一次碰撞的时刻
RV = zeros(1,len - 1);      %   每一步的速度 
for i = 2:len
    rv = randi(10000);
    rv = Av(rv);
    Tab = get_T_ab(T_ab);
    T(i) = T(i - 1) + DR(i - 1) / 1e5 / rv + Tab;
    RV(i - 1) = rv;
end

%   计算单孔时间
DT = [];    %   单孔时间
Dtp = [];   %   越孔位置对应的DR/T序号
V = [];     %   越孔速度大小
Vth = [];   %   越孔速度方向角
kxx = kx(1);
kyy = ky(1);
Dtp(1) = 1;
j = 1;
for i=2:len
    if(kx(i)~=kxx||ky(i)~=kyy)
        kxx = kx(i);
        kyy = ky(i);
        j = j + 1;
        Dtp(j) = i;
        V(j) = RV(i - 1);
        Vth(j) = TH(i - 1);
    end
end
for i = 1:j-1
    DT(i) = T(Dtp(i+1)) - T(Dtp(i));
end
%   略去开头结束两个孔
DT = DT(2:j-2);    
V = V(3:j-1);
Vth = Vth(3:j-1);
meanlenDT = (T(Dtp(j)) -  T(Dtp(2)))/(Dtp(j) - Dtp(2));%    总时间/碰撞次数=平均每一步的时间
Npore = length(DT);%    孔的数量

end
function [DR, kx, ky, TH, eX,eY] = per_walk(filename, R, a)
%   输入：路径文件，(R,a)几何参数
%   输出：DR每条路径的长度，(kx,ky)每个路径点对应的孔编号,TH每条路径的角度;(eX,eY)有效路径点
A = readtable(filename);
X = table2array(A(:,2));%  x向量
Y = table2array(A(:,3));%  y向量
%   每一步步长
DX = [];    %   每一步走的X分量
DY = [];    %   每一步走的Y分量
DR = [];    %   每一步走的长度
TH = [];    %   每一步的方向角
eX = [];    %   实际的碰撞点
eY = [];

eX(1) = X(1);
eY(1) = Y(1);
len = length(X);
j = 1;
for i = 2:len
    DX(j) = X(i) - X(i - 1);
    eX(j + 1) = X(i);   %  由于程序的原因，粒子可能在同一个位置记录两次
    DY(j) = Y(i) - Y(i - 1);
    eY(j+1) = Y(i);
    DR = (DX(j)^2+DY(j)^2)^0.5;
    if(abs(DR)<0.01)    %    排除多记录的位置
        continue;
    end
    j = j + 1;
end
DX = DX(1:j-1);
DY = DY(1:j-1);
eX = eX(1:j);
eY = eY(1:j);
DR = DX.^2 + DY.^2;
DR = DR.^0.5;
TH = atan2(DY,DX);
%   每个碰撞点的位置编号
len = length(eX);
sd = [];%   某个位置的参数
for i = 1:len
    sd(i,:) = scene_porethroat([eX(i) eY(i)], R, a);
end
kx = sd(:,2);%  表示所在孔的圆心的编号
ky = sd(:,3);
end
function s = get_T_ab(T)
%   给定平均时间T，得到一个指数分布的吸附时间
if(T < 1e-11)
    s = 0;
end
s = -T*log(rand());
end
function th0 = get_th(TH, kx, ky, kx1, ky1)
%   假设TH为跨越边界的某一步的方向角，kxky是起点所在圆标号，kx1ky1是终点所在圆标号
%   计算某一步相对于边界的角度，边界法向方向角记为0
%   th0为修正后的方向角
if(kx == kx1 && ky == ky1)              % 没有跨越边界
    th0 = 0;
    return;
end
if(abs(kx-kx1) == 2 && ky == ky1)       %  水平跨越
    th0 = TH;
elseif(abs(kx - kx1) == 1 && abs(ky - ky1) ==1)
    if((kx - kx1) - (ky - ky1) == 0)    % 斜向上
        th0 = TH - pi/3;
    else                                %  斜向下
        th0 = TH + pi/3;
    end
else                                    %   或许存在的其他情况
    th0 = 0;
end
end
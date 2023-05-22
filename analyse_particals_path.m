%   由路径库,赋予运动速度，满足2D玻尔兹曼分布
%   考虑有分布的吸附时间，时间为指数分布
%   画出均方位移与时间的图，得到扩散系数
%   记录了吸附时间-扩散系数的数据
%   记录了粒子数-扩散系数的数据，与REV有关。
clc;
clear;
filenamesta = 'data\path_six\path_2000_six_R_1.1_a_0.25\rxT_circle_onepartical_sixcir_';
filenameend = '_.txt';

NX = [];%   索引时间的位置
NY = [];
X = [];%    碰撞时刻的位置
Y = [];
DX = [];%   单次碰撞的长度
DY = [];
DR = [];
eX = [];%   每次碰撞所在位置（X有重复）
eY = [];

T = [];%    时间

sum_of_m = 1000;%    路径库文件总数
Tab = 0;%    平均吸附时间，运动特征时间为1e-8

filename_v = 'data\v_10000\Blotzmann_water_T_40.txt';
Av = readtable(filename_v);
Av = table2array(Av);

for m = 1:sum_of_m
    filenamemid = int2str(m - sum_of_m * floor((m - 1) / sum_of_m));
    filename = strcat(filenamesta,filenamemid,filenameend);

    A = readtable(filename);
    X = table2array(A(1:end,2));%  x向量
    Y = table2array(A(1:end,3));%  y向量

    j = 1;
    len = length(X);
    eX(1,m) = X(1);
    eY(1,m) = Y(1);
    for i = 2:len
        DX(j) = X(i) - X(i - 1);
        DY(j) = Y(i) - Y(i - 1);
        eX(j + 1, m) = X(i);
        eY(j + 1, m) = Y(i);
        j = j + 1;
    end
    DR = DX.^2 + DY.^2;
    DR = DR.^0.5;

    len = length(DX);

    %   赋予时间，考虑吸附
    T(1,m) = 0;%    入射时间
    Ta(1,m) = 0;%  出射时间
    for i = 2:len + 1
        rv = randi(10000);
        rv = Av(rv);
        T(i, m) = Ta(i - 1, m) + DR(i - 1) / rv / 1e5;%   一个圆看作是10微米
        T_ab = get_T_ab(Tab);
%         T_ab = 0;
        Ta(i,m) = T(i,m) + T_ab;
    end
end

%   生成索引
Tmax = min(max(Ta)) * 0.99;
num_DT = 100;%  时间的分段数量
DT(1) = 0;
for i = 2:num_DT + 1
    DT(i) = DT(1) + (i - 1) * Tmax / num_DT;
end

%   查找对应位置
NX(1,:) = eX(1,:);
NY(1,:) = eY(1,:);
for m = 1:sum_of_m
    for i = 2:num_DT + 1
        row = find(Ta(:,m) > DT(i), 1);
        if(DT(i) < T(row,m))
            f = (T(row, m) - DT(i)) / (T(row, m) - Ta(row - 1, m));
            NX(i,m) = eX(row, m) - f * (eX(row, m) - eX(row - 1, m));
            NY(i,m) = eY(row, m) - f * (eY(row, m) - eY(row - 1, m));
        else
            NX(i,m) = eX(row,m);
            NY(i,m) = eY(row,m);
        end
    end
end

%   画出均方位移和时间的关系
for i=1:num_DT + 1
    NXmean(i) = sum(NX(i,:).^2)/sum_of_m;
    NYmean(i) = sum(NY(i,:).^2)/sum_of_m;
    Nmean(i) = sum(NX(i,:).^2 + NY(i,:).^2)/sum_of_m;%    <r^2>
end
figure;
plot(DT,Nmean./1e10,'DisplayName','蒙特卡洛模拟值');
title('均方位移-时间','fontsize',14);
xlabel('时间','fontsize',14);
ylabel('均方位移','fontsize',14)
px=polyfit(DT,NXmean,1);
py=polyfit(DT,NYmean,1);
p=polyfit(DT,Nmean,1);
p=p(1) / 1e10;
px=px(1) / 1e10;
py=py(1) / 1e10;
%%
%   画出扩散系数-吸附时间图
T = [1e-11 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5];
D = [3.0619e-04 3.0651e-04 3.0239e-04 2.7330e-04 1.3709e-04 2.2893e-05 2.4522e-06];
% figure;
hold on;
plot(T,D./4);
title('扩散系数-平均吸附时间','fontsize',14);
xlabel('平均吸附时间','fontsize',14);
ylabel('扩散系数','FontSize',14);
set(gca,'xscale','log');
set(gca,'yscale','log');
box on;
plot([1e-6 5e-7],[9e-7 2e-6],'r');
plot([5e-7 5e-7],[9e-7 2e-6],'r');
plot([1e-6 5e-7],[9e-7 9e-7],'r');
patch([ss ss ss1 ss1], [6.1305e-7 1e-4 1e-4 6.1305e-7],'r','facealpha',0.2)
%%  画出粒子数-扩散系数图
x = [1 2 5 10 15 20 21 25 30 40 50 60 70 80 90 100 200 500];
y = [2.8637e-4 4.8535e-4 4.2222e-4 3.1873e-4 3.6257e-4 3.2977e-4 3.2572e-04 3.1409e-4 2.9578e-4 2.9364e-4 2.7364e-4 2.7969e-4 ...
    3.0858e-4 3.0229e-4 2.9126e-4 2.9122e-4 2.9589e-4 2.9426e-4];
hold on;
plot(x,y,'r');
set(gca,'xscale','log');
title('算出的扩散系数-粒子数','fontsize',14);
xlabel('粒子数','fontsize',14);
ylabel('扩散系数','FontSize',14);
a = y(end);
plot([1 1e3],[a*0.9 a*0.9],'b')
plot([1 1e3],[a*1.1 a*1.1],'b')


function s = get_T_ab(T)
%   给定平均时间，得到一个指数分布的吸附时间
%   T为平均吸附时间
if(T < 1e-14)
    s = 0;
else
    s = -T*log(rand());
end
end
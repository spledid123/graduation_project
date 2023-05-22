%   分析均匀排列的圆形阵列弹道数据，试图找到其与随机游走模型对应的参数
%%
% analyse_onepartial_path;
%   分析粒子的轨迹，得到孔/喉的弹道数据；分别是记录每次孔内弹道的DRmax和每次喉内弹道的DRmin;所有单孔内弹道数量Dtpore,所有单喉内弹道数量Dtth;
%%
%   分析弹道长度的数据

%   处理喉道数据，喉道数据最小值大于0,大概为喉道最窄处长度
%   喉道长度形状分布大概为指数分布曲线向右移动
%   拟合的函数是f(x)=exp((x-a)/lam)/lam;x>a

j = 1;
for i=1:length(DRmin)
    if(DRmin(i)>0.1)%   喉道最短处也是0.2；事实上，去掉这些计算误差后，最小值就是0.2
        DRmin(j)=DRmin(i);
        j=j+1;
    end
end
DRmin = DRmin(1:j-1);
DRmin_new = DRmin-min(DRmin);%  减掉最小值，方便指数分布拟合
figure;
histogram(DRmin_new,200);%   与指数分布的拟合很好
title('喉道弹道距离分布,归0版','FontSize',14);

%   处理孔内弹道长度数据
%   孔内弹道长度同样大于0，而最大值大于孔的长度，因为粒子可能越孔
%   处理数据，方便拟合，去掉大于孔道长度的数据，减掉最小值（这里指孔的最短距离）

j=1;
for i=1:length(DRmax)
    if(DRmax(i)<1.9)%   这里超过1.9指的就是越孔
        DRmax(j)=DRmax(i);
        j=j+1;
    end
end
DRmax=DRmax(1:j);
min(DRmax);
DRmax_new = DRmax - min(DRmax);
figure;
histogram(DRmax_new);%   拟合结果表示，与指数分布拟合并不成功。
title('孔内弹道距离分布,归0版','FontSize',14);

%%
%   分析单孔/喉弹道数量的数据
%   弹道数量的分布大抵是几何分布，如下办法去找参数
%   P(X=k)=(1-p)^(k-1)p;lnP = k*ln(1-p)+ln(p/(1-p))
%   先画出概率表，再做线性拟合，就找到了p
%   单孔弹道数量数据
max_Dt = max(Dtpore);
Ptp = zeros(1,max_Dt);%     Ptp(k)=n(X=k),为单孔内弹道数量为k的频数
for i=1:length(Dtpore)
    a = Dtpore(i);
    a=int16(a);
    Ptp(a) = Ptp(a) +1;
end
j = 1;
for i=1:max(Dtpore)%    避免出现频数/概率为0的情况，提前去掉后半截；可能出现n(X<=9)>0,n(X=10)=0,n(X=11)=1的情况,此处j=10
    if(Ptp(i)==0)
        j = i;
        break;
    end
end
Ptp=Ptp(1:j-1);
for i=1:j-1
    Ln_Ptp(i) = log(Ptp(i)/sum(Ptp));%     Ptp=ln(P(X=k));
end
%   画图(k,log(P(X=k)))，做线性拟合
figure;
x = 1:j-1;
histogram(Dtpore);
title('单孔弹道数量分布','FontSize',14);
p = polyfit(x,Ln_Ptp,1);
coefficient = p(1);%    斜率，用来估计概率函数参数,ln(1-p)=coefficient
p_Dtpore = 1 - exp(coefficient);%  给出参数

%  单喉弹道数量数据
max_Dt = max(Dtth);
Ptp = zeros(1,max_Dt);
for i=1:length(Dtth)
    a = Dtth(i);
    a=int16(a);
    Ptp(a) = Ptp(a) +1;
end
j = 1;
for i=1:max(Dtth)%    避免出现频数/概率为0的情况，提前去掉后半截；可能出现n(X<=9)>0,n(X=10)=0,n(X=11)=1的情况,此处j=10
    if(Ptp(i)==0)
        j = i;
        break;
    end
end
Ptp=Ptp(1:j-1);
for i=1:j-1
    Ln_Ptp(i) = log(Ptp(i)/sum(Ptp));%     Ptp=ln(P(X=k));
end
%   画图(k,log(P(X=k)))，做线性拟合
figure;
x = 1:j-1;
histogram(Dtth);
title('单喉弹道数量分布','FontSize',14);
p = polyfit(x,Ln_Ptp,1);
coefficient = p(1);%    斜率，用来估计概率函数参数,ln(1-p)=coefficient
p_Dtpore = 1 - exp(coefficient);%  给出参数
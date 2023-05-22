%   分析一下单粒子轨迹，文件记录了所有碰撞点
%   第一个点是原点
%%
%   从轨迹表格中将弹道长度提取出来
%   这个小节一定是要运行的
clc;
clear;
A = readtable('data/onepartical_long_path_six\rxT_circle_onepartical_sixcir_R_1.1_a_0.25(3).txt');
X = table2array(A(:,2));%  x向量
Y = table2array(A(:,3));%  y向量
X = X(2:end);%  去除原点
Y = Y(2:end);

%   预分配内存
len = length(X);
DX = [];    %   每一步走的X分量
DY = [];    %   每一步走的Y分量
DR = [];     %   每一步走的长度
Theta = []; %   每一步的角度

eX = [];%   实际的碰撞点
eY = [];

eX(1) = X(1);
eY(1) = Y(1);

j = 1;
for i = 2:len
    DX(j) = X(i) - X(i - 1);
    eX(j+1) = X(i);%  由于程序的原因，粒子可能在同一个位置记录两次
    DY(j) = Y(i) - Y(i - 1);
    eY(j+1) = Y(i);
    DR=(DX(j)^2+DY(j)^2)^0.5;
    if(abs(DR)<0.01)%    排除多记录的位置
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
Theta = atan2(DY, DX);
Theta = Theta ./ (2*pi).*360;%  每一步轨迹的方向
figure;
histogram(DR);
%   赋予时间
%   考虑吸附
filename_v = 'data\v_10000\Blotzmann_water_T_40.txt';
Av = readtable(filename_v);
Av = table2array(Av);
T = [];
i = 1;
T(i) = 0;
len = length(eX);
for i = 2:len
    rv = randi(10000);
    rv = Av(rv);
%     Tab = get_T_ab(1e-6);
    Tab = 0;
    T(i) = T(i - 1) + DR(i - 1) / 1e5 / rv + Tab;
end

%%  cubecir,归一化，看碰撞点在圆上的分布
%   归一化，将碰撞点归一到同一个孔内
R=[1.1 1.1];
kx = floor((eX/R(1)-1)/2);
ky = floor((eY/R(2)-1)/2);
rX = eX-2*(R(1)).*(kx+1);%   归一化的坐标
rY = eY-2*(R(2)).*(ky+1);
%   归到同一个圆里
len = length(rX);
for i = 1:len
    if(rX(i) < 0)
        rX(i) = rX(i) + 2 * R(1);
    end
    if(rY(i) < 0)
        rY(i) = rY(i) + 2 * R(2);
    end
end
%   给出相对圆的角度
Phi = atan2(rY - R(2), rX - R(1));
Phi = Phi / pi * 180;
figure;
histogram(Phi,'BinLimits',[-180,180]);%   结果是均匀分布

%%  cubecir,分析弹道长度分布
% %    分析弹道的分布规律
DXmin = [];%    喉道内的弹道
DYmin = [];
DRmin = [];
DXmax = [];%    孔内的弹道，包括越孔
DYmax = [];
DRmax = [];

len = length(DR);
j = 1;
k = 1;
for i = 1:len
    if(DR(i) < 1.1113)%     靠长度区分孔和喉
        DXmin(j) = DX(i);
        DYmin(j) = DY(i);
        DRmin(j) = (DX(i)^2+DY(i)^2)^0.5;
        j = j + 1; 
    else
        DXmax(k) = DX(i);
        DYmax(k) = DY(i);
        DRmax(k) = (DX(i)^2+DY(i)^2)^0.5;
        k = k + 1;
    end
end
figure;
histogram(DRmin);
title('喉内弹道长度分布');
figure;
histogram(DRmax);
title('孔内弹道长度分布');

%%  cubecir,统计孔/喉的弹道数量分布
Dflag = [];%    记录每一段弹道是属于哪个结构
len = length(DR);
for i = 1:len
    if(DR(i) < 1.1113)
        Dflag(i) = 0;%  表示喉
    else
        Dflag(i) = 1;%  表示孔
    end
end
Dtp = [];%    每个孔的弹道数量；
Dtth = [];%     每个喉的弹道数量
j = 0;
k = 0;

if(Dflag(1)==1)
    j = j + 1;
    Dtp(j) = 1;
else
    k =  k + 1;
    Dtth(k) = 1;
end
for i=2:len
    if(Dflag(i) == Dflag(i-1))
        if(Dflag(i)==0)
            Dtth(k) = Dtth(k) + 1;
        else
            Dtp(j) = Dtp(j) + 1;
        end
    else
        if(Dflag(i)==0)
            k = k + 1;
            Dtth(k) = 1;
        else
            j = j + 1;
            Dtp(j) = 1;
        end
    end
end
%% cubecir,统计单孔停留次数和时间
R=[1.1 1.1];
kx = floor((eX/R(1)-1)/2);
ky = floor((eY/R(2)-1)/2);
len = length(eX);
kxx = kx(1);
kyy = ky(1);
DTp = [];
Dtp(1) = 1;
DTp(1) = 0;
j = 1;
for i=2:len%    获得Dp
    if(kx(i)~=kxx||ky(i)~=kyy)
        kxx = kx(i);
        kyy = ky(i);
        j = j + 1;
        Dtp(j) = i;
        DTp(j) = T(i);
    end
end

for i = 1:j-1
    Dtp(i) = Dtp(i + 1) - Dtp(i);
    DTp(i) = DTp(i+1) - DTp(i);
end
Dtp = Dtp(1:j-1);
DTp = DTp(1:j-1);

meanDp = sum(Dtp)/length(Dtp);%   11.2759
figure;
histogram(Dtp);
Ptp = zeros(1,max(Dtp));
for i=1:length(Dtp)
    Ptp(Dtp(i)) = Ptp(Dtp(i))+1;
end
sump = length(Dtp);
for i=1:length(Ptp)
    Ptp(i) = Ptp(i)/sump;
end
%%  sixcir,碰撞点的性质，以下小节依赖该小节
%   这个模型可以直接根据位置判断是在孔/喉，而cubecir只能根据长度来
len = length(eX);
sd = [];%   某个位置的参数
for i = 1:len 
    sd(i,:) = scene_porethroat([eX(i) eY(i)],[1.1 1.1],0.25);
end
kx = sd(:,2);%  表示所在孔的圆心的编号
ky = sd(:,3);
ispore = sd(:,4);%  ispore==0表示在喉，==1表示在孔
 

%%   sixcir,单孔(孔+喉)碰撞次数分布
Dtp = [];%   在某个孔的次数
len = length(eX);
kxx = kx(1);
kyy = ky(1);
Dtp(1) = 1;
j = 1;
for i=2:len%    获得Dp
    if(kx(i)~=kxx||ky(i)~=kyy)
        kxx = kx(i);
        kyy = ky(i);
        j = j + 1;
        Dtp(j) = i;
    end
end
for i = 1:j-1
    Dtp(i) = Dtp(i + 1) - Dtp(i);
end
Dtp = Dtp(1:j-1);

meanDp = sum(Dtp)/length(Dtp);%   11.2759
figure;
histogram(Dtp);
Ptp = zeros(1,max(Dtp));
for i=1:length(Dtp)
    Ptp(Dtp(i)) = Ptp(Dtp(i))+1;
end
sump = length(Dtp);
for i=1:length(Ptp)
    Ptp(i) = Ptp(i)/sump;
end
%%  sixcir,单孔时间分布
DT = [];%   单孔时间分布
Dtp = [];
len = length(eX);
kxx = kx(1);
kyy = ky(1);
Dtp(1) = 1;
j = 1;
for i=2:len%    获得Dp
    if(kx(i)~=kxx||ky(i)~=kyy)
        kxx = kx(i);
        kyy = ky(i);
        j = j + 1;
        Dtp(j) = i;
    end
end
for i = 1:j-1
    DT(i) = T(Dtp(i+1)) - T(Dtp(i));
end
DT = DT(1:j-1);

meanDT = sum(DT)/length(DT);%   11.2759
figure;
histogram(DT);


%%   sixcir,孔/喉碰撞次数/弹道数量(看作一样的)的分布
%   求概率
j = 1;
len = length(eX);
%   喉
Dtth = [];%   喉道碰撞次数
Dtth(1) = 1;
f = 0;% 上一次是否在喉，1表示在
for i=1:len
    if(ispore(i)==0 && f==0)
        f = 1;
        Dtth(j) = 1;
    elseif(ispore(i)==0 && f==1)
        Dtth(j) = Dtth(j) + 1;
    elseif(ispore(i)==1 && f==1)
        f = 0;
        j = j + 1;
    else
    end
end

meanDt = sum(Dtth)/(length(Dtth));% 4.0733
figure;
histogram(Dtth);
title('喉内弹道数量','FontSize',14);
Ptt = zeros(1,max(Dtth));
for i=1:length(Dtth)
    Ptt(Dtth(i)) = Ptt(Dtth(i))+1;
end
sump = length(Dtth);
for i=1:length(Ptt)
    Ptt(i) = Ptt(i)/sump;
end
%   孔
Dtpore = [];
tx = [];%   所有在孔的点，对应的圆心编号
ty = [];
f = [];
j = 1;
for i=1:len-1
    if(ispore(i)==1)
        tx(j) = kx(i);
        ty(j) = ky(i);
        f(j) = 0;
        if(ispore(i+1)==0)
            f(j) = 1;%  下一步到喉
        end
        j = j+1;
    end
end
tx(j) = kx(len);
ty(j) = ky(len);
f(j) = 0;

j = 1;
Dtpore(1) = 1;
len = length(tx);
kxx = tx(1);
kyy = ty(1);
for i=2:len
    if(tx(i)~=kxx||ty(i)~=kyy||f(i-1)==1)%  上一步是喉或者另一个孔
        kxx = tx(i);
        kyy = ty(i);
        j = j + 1;
        Dtpore(j) = i;
    end
end
for i = 1:j-1
    Dtpore(i) = Dtpore(i + 1) - Dtpore(i);
end
Dtpore = Dtpore(1:j-1);

figure;
histogram(Dtpore);
title('孔内弹道数量','FontSize',14);
meanDp = sum(Dtpore)/length(Dtpore);%   10.0972
Ptpo = zeros(1,max(Dtpore));%  P(k)
for i=1:length(Dtpore)
    Ptpo(Dtpore(i)) = Ptpo(Dtpore(i))+1;
end
sump = length(Dtpore);
for i=1:length(Ptpo)
    Ptpo(i) = Ptpo(i)/sump;
end
%%   sixcir,孔/喉内单步长度
DRmax = [];%    在孔内弹道长度分布
DRmin = [];%    喉内弹道长度分布
thp=[];%    孔内弹道方向分布
tht=[];%    喉内弹道方向分布
m = 1;
n = 1;
len=length(eX);
%   区分'孔内'弹道和'喉内弹道';
for i = 1:len-1
    if(ispore(i)==1)
        thp(m) = atan2(DY(i),DX(i));
        DRmax(m)=(DX(i)^2+DY(i)^2)^0.5;
        m = m + 1;
    elseif(ispore(i)==0)
        tht(n)=atan2(DY(i),DX(i));
        DRmin(n)=(DX(i)^2+DY(i)^2)^0.5;
        n = n + 1;
    end
end
figure;
histogram(DRmax);
title('孔内弹道长度分布','FontSize',14);
figure;
histogram(DRmin);
title("喉内弹道长度分布",'fontsize',14);
%   删去多余部分
j = 1;
DRMAX = [];
for i=1:length(DRmax)
    DRMAX(j) = DRmax(i);
    if((DRmax(i) > 2) == 0)
        j = j + 1;
    end
end
j = 1;
DRMIN = [];
for i=1:length(DRmin)
    DRMIN(j) = DRmin(i);
    if((DRmin(i) > 1.02) == 0)
        j = j + 1;
    end
end

function s = get_T_ab(T)
%   给定平均时间，得到一个麦克斯韦分布的吸附时间
    s = -T*log(rand());
end
    







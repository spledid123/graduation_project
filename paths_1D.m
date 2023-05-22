%   多粒子在一维直管的运动
%%   吸附时间占主导
%   可以加上运动时间
clc;
clear;
sum_of_m = 300;%    粒子数量
num_of_alp = 50000;
fcosin = @(x) asin(2 * x -1);%  余弦定律

X = [];
T = [];
NX = [];

for m = 1:sum_of_m
X(1,m) = 0;
T(1,m) = 0;
for i = 2:num_of_alp
    alpha = fcosin(rand()); 
    X(i,m) = X(i - 1, m) + tan(alpha);
    T(i,m) = T(i - 1, m) + T_absorb();
end
end

Tmax = min(max(T)) - 0.01;
num_DT = 100;
DT(1) = 0;
for i = 2:num_DT + 1
    DT(i) = DT(i - 1) + Tmax / num_DT;
end

%   查找对应位置
NX(1,:) = X(1,:);
for m = 1:sum_of_m
for i = 2:num_DT + 1
    [row, col] = find(T(:,m) > DT(i), 1);
    NX(i,m) = X(row, m);
end
end

%   画出均方位移和时间的关系
for i=1:num_DT + 1
    NXmean(i) = sum(NX(i,:).^2)/sum_of_m;
end
plot(DT,NXmean);
%%   一维直管，没有吸附
clc;
clear;

filename_v = 'D:\Download\learn\code\mat\v_10000\Blotzmann_water_T_40.txt';
Av = readtable(filename_v);
Av = table2array(Av);

sum_of_m = 300;
num_of_alp = 50000;
fcosin = @(x) asin(2 * x -1);

X = [];
T = [];
NX = [];

for m = 1:sum_of_m
    X(1,m) = 0;
    T(1,m) = 0;
    for i = 2:num_of_alp
        alpha = fcosin(rand());
        X(i,m) = X(i - 1, m) + tan(alpha);
        %   没有吸附
        rv = randi(10000);
        rv = Av(rv);
        T(i,m) = T(i - 1, m) + 1 / cos(alpha) / rv / 1e5;
    end
end

Tmax = min(max(T)) * 0.99;
num_DT = 100;
DT(1) = 0;
for i = 2:num_DT + 1
    DT(i) = DT(i - 1) + Tmax / num_DT;
end

%   查找对应位置，没有吸附
NX(1,:) = X(1,:);
for m = 1:sum_of_m
    for i = 2:num_DT + 1
        [row, col] = find(T(:,m) > DT(i), 1);
        f = (T(row, m) - DT(i)) / (T(row, m) - T(row - 1, m));
        NX(i,m) = X(row - 1, m) + f * (X(row, m) - X(row - 1, m));
    end
end

%   画出均方位移和时间的关系
for i=1:num_DT + 1
    NXmean(i) = sum(NX(i,:).^2)/sum_of_m;
end
plot(DT,NXmean);
p = polyfit(DT,NXmean,1);
coefficient = p(1) / 1e10;
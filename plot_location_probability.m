% 画出二维概率密度图

%   同一时刻多粒子的位置表
A=readtable("data\onepartical_long_path_cir\rxT_circle_onepartical_cubecir_R_1.1.txt");
X=table2array(A(1:end,2));
Y=table2array(A(1:end,3));
%   归一化
R=[2.2 2.2];
%   为只有碰撞点的数据赋予时间
filename_v = 'data\v_10000\Blotzmann_water_T_40.txt';
Av = readtable(filename_v);
Av = table2array(Av);
j = 1;
len = length(X);
eX(1) = X(1);
eY(1) = Y(1);
for i = 2:len
    DX(j) = X(i) - X(i - 1);
    DY(j) = Y(i) - Y(i - 1);
    eX(j + 1) = X(i);
    eY(j + 1) = Y(i);
    j = j + 1;
end
DR = DX.^2 + DY.^2;
DR = DR.^0.5;
len = length(DX);
%   赋予时间
T(1) = 0;%    入射时间
Ta(1) = 0;%  出射时间
for i = 2:len + 1
    rv = randi(10000);
    rv = Av(rv);
    T(i) = Ta(i - 1) + DR(i - 1) / rv / 1e5;%   一个圆看作是10微米
%     T_ab = get_T_ab(Tab);
            T_ab = 0;
    Ta(i) = T(i) + T_ab;
end
%   生成索引
Tmax = min(max(Ta)) * 0.99;
num_DT = 5e4;%  时间的分段数量
DT(1) = 0;
for i = 2:num_DT + 1
    DT(i) = DT(1) + (i - 1) * Tmax / num_DT;
end
%   查找对应位置
NX(1) = eX(1);
NY(1) = eY(1);
for i = 2:num_DT + 1
    row = find(Ta > DT(i), 1);
    if(DT(i) < T(row))
        f = (T(row) - DT(i)) / (T(row) - Ta(row - 1));
        NX(i) = eX(row) - f * (eX(row) - eX(row - 1));
        NY(i) = eY(row) - f * (eY(row) - eY(row - 1));
    else
        NX(i) = eX(row);
        NY(i) = eY(row);
    end
end
rx = NX;
ry = NY;
%%
kx = floor((rx/R(1)-1)/2);
ky = floor((ry/R(2)-1)/2);
rx1 = rx-2*(R(1))*(kx+1);
ry1 = ry-2*(R(2))*(ky+1);
rx1 = rx1';
ry1 = ry1';

%   画出二维的概率密度图
Tlim=2.2;
dr = 0.5;
gridx1=-Tlim:dr:Tlim;
gridx2=-Tlim:dr:Tlim;
[x1,x2]=meshgrid(gridx1,gridx2);
figure;
%   连续化，f为网格的概率密度
x10=x1(:);
x20=x2(:);
xi=[x10 x20];
[f,x]=ksdensity([rx1 ry1],xi);
f=reshape(f,length(gridx1),length(gridx2));
%   离散网格，f为网格粒子数量
% f=zeros(length(gridx1),length(gridx2));
% for i=1:length(rx1)
%     a = floor((rx1(i)+Tlim)/dr)+1;
%     b = floor((ry1(i)+Tlim)/dr)+1;
%     f(a,b)=f(a,b)+1;
% end
    
pcolor(x1, x2, f); 
% shading interp
set(gca,'xlim',[-Tlim,Tlim]);
set(gca,'ylim',[-Tlim,Tlim]);
colorbar;
grid off;
hold on;
% 画圆
plot_cir(-1.1,1.1,1,1);
plot_cir(-1.1,-1.1,1,1);
plot_cir(1.1,1.1,1,1);
plot_cir(1.1,-1.1,1,1);
axis equal;
xlim([-Tlim Tlim]);
ylim([-Tlim Tlim]);

%   输入：多粒子在不同时刻的位置，可能在周期性边界内
%   判断是否达到稳定
%   分析速度分布

%%  分析两个区块的粒子数量比以及数密度的比随时间的变化
% 两个区块的边界在x = 0
clear;
f = [];         % 不同时刻右边圆阵列/左边圆阵列，粒子数量比例
numr = 110;      % 右边圆阵列圆数量
numl = 500;      % 左边圆阵列圆数量
Rr = 1;         % 右边圆阵列圆半径
Rl = 1;         % 左边圆阵列圆半径
tr = 5;      % 右边圆阵列圆半径+喉道宽度一半
tl = 1.1;      % 左边圆阵列圆半径+喉道宽度一半
tsdr = 0;      %   区块与圆阵列的间隔
tsdl = 0;
%   虽然开启并行区需要不少时间，但并行确实节省时间
for i = 10:10:900
    filenamesta = 'data\dif_bulk_pore_cir_R_1_5_110_1_1.1_500_N_50000\rxT_circle_T_';
    filenamemid = num2str(1*i);
    filenameend = '_.txt';
    filename = strcat(filenamesta,filenamemid,filenameend);
    A = readtable(filename);
    len = size(A);
    len = len(1);
    %   粒子数
    path_start = 1;
    path_end = len;
    rx=table2array(A(path_start:path_end,1));
    ry=table2array(A(path_start:path_end,2));
    %   减去多孔介质(圆阵列)与区块间隔部分
    fxr= (rx >= tsdr) & (rx <= (tsdr + 2 * numr * tr)) & (ry >= tsdr) & (ry <= (tsdr + 2 * numr * tr));
    fxl= (rx <= -tsdl) & (rx >= -(tsdl + 2 * numl * tl)) & (ry >= tsdl) & (ry <= (tsdl + 2 * numl * tl));
    f(i) = sum(fxr)/sum(fxl);
end
%   左边区块是空的
%   右边区块孔隙面积/左边区块面积
% numx = 20;
% numy = 20;
% f_s = (4*numx*numy*R(1)*R(2)-numx*numy*pi)/(8*numx*numy*R(1)*R(2)-numx*numy*pi);

%   左右都是圆阵列
Sr = (numr * numr)*(4 * tr * tr - pi * Rr * Rr);%  右边孔隙面积
Sl = (numl * numl)*(4 * tl * tl - pi * Rl * Rl);%   左边孔隙面积
f_s = Sr / Sl;

figure;
hold on;
box on;
title('右边多孔介质的粒子数与左边的比例','fontsize',14);
xlabel('时间','FontSize',14);
plot(10:10:840,f(10:10:840),'DisplayName','粒子在右边多孔介质中的数量/在左边的数量');
plot([1 1000],[f_s f_s],'DisplayName','右边孔隙面积/左边孔隙面积');

%   左右粒子数密度的比
figure;
hold on;
box on;
title('分子数密度之比：右边比左边','fontsize',14);
xlabel('时间','FontSize',14);
set(gca,'ylim',[0 1.2]);
fn = f/f_s;
plot(10:10:840,fn(10:10:840),'DisplayName','右边孔隙数密度/左边孔隙数密度');
%   meanfn求平均取值范围
% fn_sta = 800;
% fn_end = 840;
% meanfn = sum(fn(fn_sta:fn_end))/(fn_end - fn_sta + 1);%   平均值：右边数密度/左边数密度
% plot([1 2000],[meanfn meanfn],'DisplayName','右边孔隙数密度/左边孔隙数密度');
%% 根据圆内多孔介质模型的特点，单纯生成随时间变化的圆内粒子/多孔介质内粒子数，以判断是否平衡
%   计算圆阵列内部以及与外部的数密度
clear;
tic;
r = 500;    %   圆的半径
Rp = 1.5;  %  阵列孔喉比
nump = 50;   % 阵列行数
R = [r Rp nump];
num_cir = zeros(5500,1);%   圆内粒子数
num_pore = zeros(5500,1);%  多孔介质内粒子数
num_1 = zeros(5500,1);
num_2 = zeros(5500,1);
num_3 = zeros(5500,1);
num_4 = zeros(5500,1);
num_5 = zeros(5500,1);
num_6 = zeros(5500,1);
flag = 0;    %   1表征这是存在多孔介质，0表示没有
% 循环遍历每个表格
ii = 10:10:500;
N_T = 1;
j = 1;
V = [];
TH = [];
tic;
for i = ii
    % 从表格中提取粒子的x和y坐标,可能有两个文件
    filenamesta = 'data\bulk_pore_pm18_R_500_N_1000000_dT_1_2000\rxT_circle_T_';
    %     filenamesta1 = 'data\bulk_pore_cir_R_500_1.1_100_N_1000000_dT_1_2000(2)\rxT_circle_T_';
    filenamemid = num2str(N_T*i);
    filenameend = '_.txt';
    filename = strcat(filenamesta,filenamemid,filenameend);
    %     filename1 = strcat(filenamesta1,filenamemid,filenameend);
    A = readtable(filename);
    %     A1 = readtable(filename1);
    %     x = [A.rx;A1.rx];
    %     y = [A.ry;A1.ry];
    %     sd = [A.sd;A1.sd];
    x = [A.rx];
    y = [A.ry];
    sd = [A.sd];
    v = A.v;
    th = atan2(A.gy,A.gx);

    %   根据圆-多孔介质模型的特点进行计数筛选
%         t_pore = (x.^2 + y.^2).^0.5 < 200;    %   半径100的圆弧
%         t_pore = (x < Rp * nump) & (x > -Rp * nump) & (y > - Rp * nump) & (y < Rp * nump);    %   阵列
%         t_pore = (x < 180) & (x > -180) & (y > - 180) & (y < 180);  % 四个正方形或者阵列或者四个圆
    %   t_pore = (x < 150) & (x > -150) & (y > - 150) & (y < 150);
%       t_pore = (x < 110) & (x > -110) & (y > - 110) & (y < 110);
%       t_pore = (x < 205) & (x > -205) & (y > - 205) & (y < 205);
%       三角
        t_pore = [];
        for k = 1:length(v)
            t_pore(k) = triangleSDF([x(k) y(k)], [-100,-100/3^0.5], [100,-100/3^0.5], [0, 200/3^0.5]) < 0;
        end
% t_pore = zeros(1e6, 1);
% for k = 1:length(v)
%     t_pore(k) = HexagonSDF([x(k) y(k)], [0 0], 200) < 0;
% end
%         t_1 = atan2(y,x) < pi/3 & atan2(y,x) > 0 & t_pore;
%         t_2 = atan2(y,x) > pi/3 & atan2(y,x) < 2*pi/3 & t_pore;
%         t_3 = atan2(y,x) > 2*pi/3 & atan2(y,x) < pi & t_pore;
%         t_4 = atan2(y,x) > -pi & atan2(y,x) < -2*pi/3 & t_pore;
%         t_5 = atan2(y,x) > -2*pi/3 & atan2(y,x) < -pi/3 & t_pore;
%         t_6 = atan2(y,x) > -pi/3 & atan2(y,x) < 0 & t_pore;
    %   两圆
    %     t_pore = ((x + 150).^2 + (y).^2) < 205^2 | ((x - 150).^2 + (y).^2) < 205^2;   %   l
    %   三凹圆弧
%     r0 = [0 -150];
%     r1 = [150/2*3^0.5, 150/2];
%     x1 = [r0(1) + 125*cos(40/180*pi) r0(2) + 125*sin(40/180*pi)];
%     x2 = [r1(1) + 125*cos(2*pi/3+pi/2+50/180*pi) r1(2) + 125*sin(2*pi/3+pi/2+50/180*pi)];
%     x3 = [r1(1) + 125*cos(2*pi/3+pi/2-50/180*pi) r1(2) + 125*sin(2*pi/3+pi/2-50/180*pi)];
%     %   找"三角形"的三个顶点
%     X1 = [0 x2(2) + (-x2(1)) * (x2(2) - x1(2)) / (x2(1) - x1(1))];
%     X2 = [(x3(2) - x2(2)) * (x2(1) - x1(1)) / (x2(2) - x1(2)) + x2(1) x3(2)];
%     X3 = [-X2(1) X2(2)];
%     t_pore = zeros(1e6,1);
%     for k = 1:length(v)
%         t_pore(k) = triangleSDF([x(k) y(k)], X1, X2, X3) < 0;
%     end
%     t_pore = t_pore | (x.^2 + (y + 150).^2 < 125^2) | ((x-r1(1)).^2 + (y-r1(2)).^2 < 125^2)...
%                     | ((x+r1(1)).^2 + (y-r1(2)).^2 < 125^2);
%   三凸圆弧
%     t_pore = zeros(1e6,1);
%     for k = 1:length(v)
%         t_pore(k) = triangleSDF([x(k) y(k)], [-100*3^0.5 100], [100*3^0.5 100], [0 -200]) < 0;
%     end
%     t_pore = t_pore & ((x).^2 + (y + 200).^2 > 155^2) & ((x - 100*3^0.5).^2 + (y - 100).^2 > 155^2)...
%                     & ((x + 100*3^0.5).^2 + (y - 100).^2 > 155^2);



    num_pore(i) = sum(t_pore);
    t_cir = (x.^2 + y.^2).^0.5 < r-1e-3;
    num_cir(i) = sum(t_cir);
%     num_1(i) = sum(t_1);
%     num_2(i) = sum(t_2);
%     num_3(i) = sum(t_3);
%     num_4(i) = sum(t_3);
%     num_5(i) = sum(t_3);
%     num_6(i) = sum(t_3);
%     if(mod(i,10) == 0)  
        fprintf("%d\n",i);
        for k = 1:length(v)%    检验多孔介质内速度分布
            if(t_pore(k) == 1)
                V(j) = v(k);
                TH(j) = th(k);
                j = j + 1;
            end
        end
%     end
end
%   统计圆内/圆阵列内粒子速度，看分布，当然这时候上面不应该循环
% v = A.v;
% j = 1;
% for i = 1:length(v)
%     if(t_cir(i) == 1)
%         V(j) = v(i);
%         j=j+1;
%     end
% end
toc;
%  多孔介质内孔的面积 
%   阵列
% if(flag == 0)
%     Sp = (2 * nump * Rp)^2;
% else
%     Sp = (2 * nump * Rp)^2 - nump * nump * pi;
% end
% Sp = pi*(100)^2;%   半径100圆弧
%  随机圆
% Sp = (205*2)^2 * 0.7584;
% Sp = 420 ^ 2 - 4 *pi*(100)^2;
% Sp = 300^2 - 4 *50^2;
% Sp = 115*3^0.5*115*3 - 3 * 84 * 3^0.5 * 10 * 4;       %   三角
%  两圆
% th1 = acos(150/205);
% Sp = 2 * (pi * 205^2 - th1 * 205^2) + 2 * 150 * sin(th1) * 205 - 2 * pi * (205^2 - 200^2) * (130/180) - pi * 2.5^2;
%    Sp = (220)^2 - pi * (100)^2;
%   三凹圆
% Sp = 260/360*3*pi*(125)^2 + 3 * Sf(r0, x1, [-x1(1) x1(2)]) + Sf(X1, X2, X3)...
%     - 3*Sf(X1,x1,[-x1(1) x1(2)]) - 3*(260/360)*pi*(125^2-120^2);
%   三凸圆
% Sp = Sf([-100*3^0.5 100], [100*3^0.5 100], [0 -200]) - 1/2*pi*(155)^2;
% Sp = 360^2 - 4*30^2;
% Sp = (2 * nump * Rp)^2 - nump*nump*4;%`  正方形阵列
% Sp = 40011;%pm1
%   pm2
% th = 5/180*pi;
% a = 200 * (sin(th) * 3^0.5 + cos(th));
% S1 = 3^0.5*(200*sin(th))^2;
% Shu = th*200^2;
% Ssan = 200^2*sin(th)*cos(th);
% S2 = Shu - Ssan;
% S3 = S1 - S2;
% S4 = (3*a/2)*(3^0.5*a/2);
% Sp = S4 - 3*S3;
% Sp = 67059;%与精确值有0.3%的差距
%   pm3
% Sp = 41440;
%pm4
% Sp = 42434;
%pm5
% Sp = 47023.7;
% pm6
% Sp = 86165;
% pm7
% Sp = 84870;
% pm8
% Sp = 94047;
%   pm10
% Sp = 13813;
%   pm11
% Sp = 14145;
%   pm12
% Sp = 15675;
%   pm13
% Sp = 14361;
%   pm14
% Sp = 13813;
%   pm15
% Sp = 16238;
% pm16
% Sp = 82878.5;
% %   pm17
% Sp = 15674.7;
% %   pm18
Sp = 16612;



f_p = num_pore / Sp; %  多孔介质数密度
%   计算方阵外数密度
% Sc = pi * (r - 1e-3) ^ 2 - (2 * nump * Rp)^2;
% Sc = pi * (r - 1e-3)^2 - (360)^2;
% Sc = pi * (r - 1e-3) ^ 2 -115*3^0.5*115*3;
% Sc = pi * (r - 1e-3) ^ 2 - (2 * (pi * 205^2 - th1 * 205^2) + 2 * 150 * sin(th1) * 205);
% Sc = pi * (r - 1e-3)^2 - (220)^2 - 3*pi*100^2;
%   三凹圆
%   Sc = pi * (r - 1e-3)^2 - (260/360*3*pi*(125)^2 + 3 * Sf(r0, x1, [-x1(1) x1(2)]) + Sf(X1, X2, X3)...
%                          - 3*Sf(X1,x1,[-x1(1) x1(2)]));
%   三凸圆
% Sc = pi * (r - 1e-3)^2 - Sf([-100*3^0.5 100], [100*3^0.5 100], [0 -200]) + 1/2*pi*(150)^2;
% Sc = pi * (r - 1e-3) ^ 2 - (2 * 205)^2;
% Sc = pi* (r-1e-3)^2 - pi*200^2;
% Sc = pi* (r-1e-3)^2 - 100*3^0.5*300;
Sc = pi * (r-1e-3)^2 - 100^2*3^0.5;
% Sc = pi * (r-1e-3)^2 - 3^1.5/2*200^2;

f_c = (num_cir - num_pore) / Sc; %  多孔介质外数密度
%  画图
figure;
hold on;
box on;
title('数密度-时间','fontsize',14);
xlabel('时间','fontsize',14);
ylabel('数密度','fontsize',14);
plot((ii).*N_T,f_p(ii),'displayname','多孔介质内数密度');
plot((ii).*N_T,f_c(ii),'displayname','多孔介质外数密度');
figure;
box on;
hold on;
title('粒子数','FontSize',14);
xlabel('时间','FontSize',14);
ylabel('粒子数','FontSize',14);
plot((ii).*N_T,num_pore(ii),'displayname','多孔介质内粒子数');
plot((ii).*N_T,num_cir(ii) - num_pore(ii),'displayname','多孔介质外数密度');
% figure;
% hold on;
% title('粒子数','FontSize',14);
% xlabel('时间','FontSize',14);
% ylabel('粒子数','FontSize',14);
% box on;
% plot(ii.*N_T, num_1(ii), 'DisplayName','1');
% plot(ii.*N_T, num_2(ii), 'DisplayName','2');
% plot(ii.*N_T, num_3(ii), 'DisplayName','3');
% plot(ii.*N_T, num_4(ii), 'DisplayName','4');
% plot(ii.*N_T, num_5(ii), 'DisplayName','5');
% plot(ii.*N_T, num_6(ii), 'DisplayName','6');
mailme([],'down',['花费了' num2str(toc/3600) '小时']);
%%  根据粒子的图表，形成速度分布
filename = 'data\bulk_pore_cir_R_500_N_20000_dT_1_200\rxT_circle_T_100_.txt';
A = readtable(filename);
x = A.rx;
y = A.ry;
v = A.v;
vx = A.gx;
vy = A.gy;
TH = atan2(vy,vx);

vv = [];
TTH = [];
j = 1;
R = 500;
%   根据圆-多孔介质模型的特点进行计数筛选
t = (x.^2 + y.^2 - R ^ 2) < -1e-3;
num_cir = sum(t);
for i = 1:length(x)
    if(t(i) == 1)
        vv(j) = v(i);
        TTH(j) = TH(i);
        j = j + 1;
    end
end
%%  对左右区块纵向分块统计，粒子的速度/数密度分布，包括方向角与速度大小的分布，按x方向将区块分割为数个部分分别统计
clear;
filename = 'data\dif_bulk_pore_cir_R_2_2.2_50_1_2.2_50_N_10000\rxT_circle_T_100_.txt';
%   x = 0左右两边的区块,圆阵列,如果一边是空白，numl = 0,如果区块的左下角(右下角)不是(0,0),调整tsdr/tsdl
numr = 50;      % 右边圆阵列圆数量
numl = 50;      % 左边圆阵列圆数量
Rr = 2;         % 右边圆阵列圆半径
Rl = 1;         % 左边圆阵列圆半径
tr = 2.2;      % 右边圆阵列圆半径+喉道宽度一半
tl = 2.2;      % 左边圆阵列圆半径+喉道宽度一半
tsdr = 0;      %   区块与圆阵列的间隔
tsdl = 0;
A = readtable(filename);
len = size(A);
len = len(1);
%   粒子数
path_start = 1;
path_end = len;
rx = table2array(A(path_start:path_end,1));
ry = table2array(A(path_start:path_end,2));
vx = table2array(A(path_start:path_end,3));
vy = table2array(A(path_start:path_end,4));
v = table2array(A(path_start:path_end,5));
TH = atan2(vy,vx);

th = v; %   进行统计的对象，可以是速度v，速度方向角TH;
%   按x方向分为数个部分,左右分开
%   如果不需要分块，则len_right = 2 * tr * numr;len_left = 2 * tl * numl;
len_right = 2 * tr * numr; %    按x方向分割区块，右边每部分的宽度
len_left = 2 * tl * numl;  %   分割区块，左边每部分的宽度
num_right_gap = zeros(1, round(2 * tr * numr / len_right)); %   右边多孔介质纵向分块粒子数
num_left_gap = zeros(1, round(2 *tl * numl / len_left));  %  左边多孔介质纵向分块粒子数
th_right_gap = cell(1, round(2 * tr * numr / len_right));    %   右边多孔介质纵向分块的统计量
th_left_gap = cell(1, round(2 * tl * numl / len_left));      %     左边多孔介质纵向分块的统计量
num_right = 0; %   右边多孔介质粒子数
num_left = 0;  %   左边多孔介质粒子数
th_right = [];%   右边多孔介质的统计量
th_left = [];%   左边多孔介质的统计量
for i = 1:len
    if((ry(i) > tsdr) && (ry(i) < tsdr + 2 * tr * numr) && (rx(i) > tsdr) && (rx(i) < tsdr + 2 * tr * numr))
        m = floor((rx(i) - tsdr)/len_right) + 1;
        num_right_gap(m) = num_right_gap(m) + 1;
        th_right_gap{1, m}(num_right_gap(m)) = th(i);
    elseif((ry(i) > tsdl) && (ry(i) < tsdl + 2 * tl * numl) && (rx(i) < -tsdl) && (rx(i) > -(tsdl + 2 * tl * numl)))
        m = floor(abs((rx(i)+tsdl))/len_left) + 1;
        num_left_gap(m) = num_left_gap(m) + 1;
        th_left_gap{1, m}(num_left_gap(m)) = th(i);
    end
end
num_right = sum(num_right_gap);
num_left = sum(num_left_gap);
th_right = table2array(cell2table(th_right_gap));
th_left = table2array(cell2table(th_left_gap));
%   数密度与x的关系
figure;
box on;
hold on;
title('纵向分块数密度','FontSize',14);
xlabel('x','fontsize',14);
ylabel('数密度','fontsize',14);
Sr = (numr * numr)*(4 * tr * tr - pi * Rr * Rr);%  右边孔隙面积
Sl = (numl * numl)*(4 * tl * tl - pi * Rl * Rl);%   左边孔隙面积
plot((1:round(2 * tr * numr / len_right)).*len_right, num_right_gap./(Sr / round(2 * tr * numr / len_right)),'o--','DisplayName','右边多孔介质纵向分块数密度');
plot((-1:-1:-round(2 * tl * numl / len_left)).*len_left, num_left_gap./(Sl / round(2 * tl * numl / len_left)),'o--','DisplayName','左边多孔介质纵向分块数密度');

%   对某个区块的速度进行统计
th = th_right;%    统计对象，可以是v,也可以是某个部分的v
[xx, yy] = plot_distribution(th, 0, 100);
figure;
box on;
hold on;
title('速度大小分布','fontsize',14);
xlabel('速度大小','fontsize',14);
ylabel('概率密度','FontSize',14);
plot(xx,yy,'displayname','速度大小分布');
%   画出速度的玻尔兹曼分布
k = 1.380649e-23;
NA = 6.02e23;
m = 18e-3/NA;
T = 40;
A = m/2/k/T;
f = @(x)2.*A.*x.*exp(-A.*x.^2);
x1 = 0:1e-2:670;
x2 = 671:1e3;
x = [x1 x2];
y = zeros(1,length(x));
y = f(x);
plot(x,y,'displayname','2D玻尔兹曼分布');

%   对某个区块的方向角进行统计
th = th_left;%   统计对象，可以是TH或者某个部分的TH
figure;
polaraxes;
hold on;
title('方向角分布','FontSize',14);
polarhistogram(th, 100,'Normalization','pdf');



function s = Sf(X1, X2, X3)
%   根据三角形三个顶点的坐标计算面积
a = norm(X2 - X1);
b = norm(X3 - X2);
c = norm(X2 - X3);
p = (a + b + c) / 2;
s = (p * (p - a) * (p - b) * (p - c))^0.5;
end
function s = Sseg(r, th)
%   给定半径与角度，计算弓形面积
s = th/2*r^2 - r^2*sin(th/2)*cos(th/2);
end





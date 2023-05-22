%   输入：多粒子在不同时刻的位置，可能在周期性边界内
%   整个周期性边界大致为九宫格：存在两种区块，分别有不同的形态，任意一种的上下左右都是另一种
%   多孔介质由圆阵列代表
%%  画粒子
filename = 'data\bulk_pore_box_R_500_2_100_N_1000000_dT_1_2000\rxT_circle_T_200_.txt';
A = readtable(filename);
len = size(A);
len = len(1);
%   粒子数
path_start = 1;
path_end = len;
rx=table2array(A(path_start:path_end,1));
ry=table2array(A(path_start:path_end,2));
figure;
plot(rx,ry,'b.');
axis equal;
hold on;
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


%% 根据一系列粒子位置的表格文件，生成粒子运动的gif,保留红色的圆形
%   该程序的逻辑是：gif的生成来自frame生成的img,而我们可以从画布gcf(也就是程序运行时的显示)生成frame.
%   我这里的frame其实是两个画布合并而来，做了一些处理——也就是说，显示的画布和gif其实不同。
%   同时生成随时间变化的圆内粒子/多孔介质内粒子数，以判断是否平衡
clear;
FileName = 'data\pic\space.gif';% 创建文件名
axis off;
axis equal;
hold on;
plot_cir_1(0,0,500);
% plot([200 200 -200 -200 200],[200 -200 -200 200 200],'r--','LineWidth',2);
%   x = 0左右两边的区块,圆阵列,如果一边是空白，numl = 0,如果区块的左下角(右下角)不是(0,0),调整tsdr/tsdl
% numr = 100;    % 右边圆阵列圆数量
% numl = 0;    % 左边圆阵列圆数量
% Rr = 1;      % 右边圆阵列圆半径
% Rl = 1;      % 左边圆阵列圆半径
% tr = 2;      % 右边圆阵列圆半径+喉道宽度一半
% tl = 2;      % 左边圆阵列圆半径+喉道宽度一半
% tsdr = -200;      %   区块间的间隔
% tsdl = 0;
% kx = 0:numr - 1;%   右边
% ky = 0:numr - 1;
% for i=1:numr
%     x0(i)=tr*(2*kx(i)+1) + tsdr;
% end
% for i=1:numr
%     y0(i)=tr*(2*ky(i)+1) + tsdr;
% end
% for i=1:numr
%     for j=1:numr
%         plot_cir(x0(i),y0(j),Rr);
%     end
% end
% kx = 0:numl - 1;%   左边
% ky = 0:numl - 1;
% for i=1:numl
%     x0(i) = - (tl*(2*kx(i)+1)) - tsdl;
% end
% for i=1:numl
%     y0(i) = tl*(2*ky(i)+1) + tsdl;
% end
% for i=1:numl
%     for j=1:numl
%         plot_cir(x0(i),y0(j),Rl);
%     end
% end
R = 500;    %   圆的半径
xlim([-R,R]);
ylim([-R,R]);
% c = load('data\cir.mat');
% c = c.cir;
% for i = 1:1300
%     plot_cir(c(i,1),c(i,2),c(i,3));
% end
% plot([-106 * 3^0.5, 106 * 3^0.5],[-106 -106],'r');
% plot([-106 * 3^0.5, 0],[-106 212],'r');
% plot([0, 106 * 3^0.5],[212 -105],'r');
frame_cir = getframe(gcf);
cla();% 清除画布(gcf)图像
num_cir = zeros(1000,1);
num_pore = zeros(1000,1);
N_T = 1;
% 循环遍历每个表格
for i = 1:1:30
    % 从表格中提取粒子的x和y坐标
    filenamesta = 'data\bulk_R_500_N_200000_dT_2000\rxT_circle_T_';
    filenamemid = num2str(N_T*i);
    filenameend = '_.txt';
    filename = strcat(filenamesta,filenamemid,filenameend);
    A = readtable(filename);
    x = A.rx;
    y = A.ry;
    %   根据圆-多孔介质模型的特点进行计数筛选
%     t = (x < tr * numr) & (x > -tr * numr) & (y < tr * numr) & (y > -tr * numr);
%     num_pore(i) = sum(t);
    t = (x.^2 + y.^2 - R ^ 2) < -1e-3;
    num_cir(i) = sum(t);
    x = x .* t;
    y = y .* t;
    % 绘制粒子图
    plot(x,y,'b.');
    plot(0,0,'w.');
    box on;
    title(['Time = ',num2str(i*N_T)])
    % 获取当前帧
    frame = getframe(gcf);
    % 将帧转换为图像
        frame.cdata = min(frame.cdata,frame_cir.cdata);%    保留圆的像素
    img = frame2im(frame);
    % 将图像转换为索引颜色
    [imgind,cm] = rgb2ind(img,256);
    cla();% 清除画布图像
    % 写入GIF文件
    if i == 1
        % 对于第一帧，开始“LoopCount”
        imwrite(imgind,cm,FileName,'gif','LoopCount',Inf,'DelayTime',0.01);
    else
        % 对于其他帧，追加到文件中
        imwrite(imgind,cm,FileName,'gif','WriteMode','append','DelayTime',0.01);
    end
end
%% 根据圆内多孔介质模型的特点，单纯生成随时间变化的圆内粒子/多孔介质内粒子数，以判断是否平衡
%   计算圆阵列内部以及与外部的数密度
clear;
r = 500;    %   圆的半径
Rp = 1.5;  %  阵列孔喉比
nump = 50;   % 阵列行数
R = [r Rp nump];
num_cir = zeros(5500,1);%   圆内粒子数
num_pore = zeros(5500,1);%  多孔介质内粒子数
flag = 0;    %   1表征这是存在多孔介质，0表示没有
% 循环遍历每个表格
ii = 10:10:370;
N_T = 1;
j = 1;
V = [];
TH = [];
tic;
for i = ii
    % 从表格中提取粒子的x和y坐标,可能有两个文件
    filenamesta = 'data\bulk_pore_box_R_500_1.5_50_N_1000000_dT_1_2000\rxT_circle_T_';
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
%         t_pore = (x.^2 + y.^2).^0.5 < 100;    %   半径100的圆弧
        t_pore = (x < Rp * nump) & (x > -Rp * nump) & (y > - Rp * nump) & (y < Rp * nump);    %   阵列
%         t_pore = (x < 180) & (x > -180) & (y > - 180) & (y < 180);  % 四个正方形或者阵列或者四个圆
    %   t_pore = (x < 150) & (x > -150) & (y > - 150) & (y < 150);
    %   t_pore = (x < 110) & (x > -110) & (y > - 110) & (y < 110);
    %   三角
    %     t_pore = [];
    %     for k = 1:length(v)
    %         t_pore(k) = triangleSDF([x(k) y(k)], [-115*3^0.5,-115], [115*3^0.5,-115], [0, 230]) < 0;
    %     end
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
% Sp = (2 * nump * Rp)^2 * 0.7584;
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
Sp = (2 * nump * Rp)^2 - nump*nump*4;%`  正方形阵列


f_p = num_pore / Sp; %  多孔介质数密度
%   计算方阵外数密度
Sc = pi * (r - 1e-3) ^ 2 - (2 * nump * Rp)^2;
% Sc = pi * (r - 1e-3)^2 - (360)^2;
% Sc = pi * (r - 1e-3) ^ 2 -115*3^0.5*115*3;
% Sc = pi * (r - 1e-3) ^ 2 - (2 * (pi * 205^2 - th1 * 205^2) + 2 * 150 * sin(th1) * 205);
% Sc = pi * (r - 1e-3)^2 - (220)^2 - 3*pi*100^2;
%   三凹圆
%   Sc = pi * (r - 1e-3)^2 - (260/360*3*pi*(125)^2 + 3 * Sf(r0, x1, [-x1(1) x1(2)]) + Sf(X1, X2, X3)...
%                          - 3*Sf(X1,x1,[-x1(1) x1(2)]));
%   三凸圆
% Sc = pi * (r - 1e-3)^2 - Sf([-100*3^0.5 100], [100*3^0.5 100], [0 -200]) + 1/2*pi*(150)^2;


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


%%    实现速度云图
filenamesta = 'data\diffusion_R_1e-5_1.1_dT_1e-5\rxT_circle_T_';
filenamemid = num2str(0.00041);
filenameend = '_.txt';
filename = strcat(filenamesta,filenamemid,filenameend);
A = readtable(filename);
r = 500;
Rp = 1.1;
nump = 30;
x = A.rx; % 提取x坐标
y = A.ry; % 提取y坐标
x = x*1e5;y = y*1e5;
kx = floor((x/Rp-1)/2);
ky = floor((y/Rp-1)/2);
X = x-2*(Rp)*(kx+1);
Y = y-2*(Rp)*(ky+1);
V = A.v; % 提取速度大小
% t_pore = (x < Rp * nump) & (x > -Rp * nump) & (y < Rp * nump) & (y > -Rp * nump);
% % t_pore = [];
% % for k = 1:length(v)
% %     t_pore(k) = triangleSDF([x(k) y(k)], [-115*3^0.5,-115], [115*3^0.5,-115], [0, 230]) < 0;
% % end
% t_cir = (x.^2 + y.^2).^0.5 < 60-1e-3;
% j = 1;
% for i = 1:length(v)
%     if(t_cir(i) == 1)
%         V(j) = v(i);
%         X(j) = x(i);
%         Y(j) = y(i);
%         j=j+1;
%     end
% end
%归一化
% t = 2 * 1.5;
% tx = t;
% ty = t;
% kx = floor(x/tx);
% ky = floor(y/ty);
% x = x - kx * tx;
% y = y - ky * ty;
% 计算不同位置的粒子平均速度
[x,y] = meshgrid(min(X):0.02:max(X),min(Y):0.02:max(Y)); % 生成方格内的网格点
v = griddata(X,Y,V,x,y); % 使用griddata函数插值得到每个网格点的速度大小
% 绘制云图，使用颜色表示速度大小
pcolor(x,y,v); % 使用pcolor函数绘制云图
shading interp; % 平滑颜色过渡
colorbar; % 显示颜色条
axis equal;
hold on;
xlabel('x'); % 标注x轴
ylabel('y'); % 标注y轴
title('Average velocity of particles at different positions'); % 标注标题
Tlim = 1.1;
plot_cir(Tlim,Tlim,1);
plot_cir(Tlim,-Tlim,1);
plot_cir(-Tlim,Tlim,1);
plot_cir(-Tlim,-Tlim,1);
xlim([-Tlim Tlim]);
ylim([-Tlim Tlim]);
% plot_cir(0,0,500);
% numr = 30;    % 右边圆阵列圆数量
% numl = 50;    % 左边圆阵列圆数量
% Rr = 1;       % 右边圆阵列圆半径
% Rl = 1.5;      % 左边圆阵列圆半径
% tr = 2;      % 右边圆阵列圆半径+喉道宽度一半
% tl = 2.2;      % 左边圆阵列圆半径+喉道宽度一半
% tsdr = -60;      %  圆阵列所在正方形与区块的间隔
% tsdl = 0;
% %   右边区块
% kx = 0:numr - 1;
% ky = 0:numr - 1;
% for i=1:numr
%     x0(i)=tr*(2*kx(i)+1) + tsdr;
% end
% for i=1:numr
%     y0(i)=tr*(2*ky(i)+1) + tsdr;
% end
% for i=1:numr
%     for j=1:numr
%         plot_cir(x0(i),y0(j),Rr);
%     end
% end
% xlim([-500 500]);
% ylim([-500 500]);
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
%%  根据多例子轨迹，画<R^2>-t，计算扩散系数
%   D = <R^2>/t
%   线性回归
filenamesta = "data\diffusion_R_1e-5_1.1_dT_1e-5\rxT_circle_T_";
filenameend = "_.txt";
N_T = 1e-5;
ii = 1:43;
r2x = zeros(1,1e3);
r2y = zeros(1,1e3);
r2R = zeros(1,1e3);
th = pi / 4; %    计算某个特定方向角的均方位移
r2thx = zeros(1,1e3);
r2thy = zeros(1,1e3);
tic;
for i = ii
    filename = strcat(filenamesta,num2str(i * N_T),filenameend);
    A = readtable(filename);
    x = A.rx;
    y = A.ry;
    len = length(x);
    r2x(i) = sum(x.^2) / len;% x方向均方位移
    r2y(i) = sum(y.^2) / len;% y方向均方位移
    r2R(i) = r2x(i) + r2y(i);      % 均方位移
    A = [x y] * [cos(th) -sin(th);sin(th) cos(th)];
    x = A(:,1);
    y = A(:,2);
    r2thx(i) = sum(x.^2) / len;%   th方向角的均方位移
    r2thy(i) = sum(y.^2) / len;%   th方向角垂直方向的均方位移
end
toc;
figure;
hold on;
box on;
xlabel('时间/s','fontsize',14);
ylabel('均方位移/m^2','fontsize',14);
title('均方位移-时间','fontsize',14);
plot(ii*N_T,r2R(ii),'DisplayName','均方位移');
plot(ii*N_T,r2x(ii),'DisplayName','x方向');
plot(ii*N_T,r2y(ii),'DisplayName','y方向');
plot(ii*N_T,r2thx(ii),'DisplayName','th方向');
plot(ii*N_T,r2thy(ii),'DisplayName','th垂直方向');
px = polyfit(ii*N_T,r2x(ii),1);
py = polyfit(ii*N_T,r2y(ii),1);
pthx = polyfit(ii*N_T,r2thx(ii),1);
pthy = polyfit(ii*N_T,r2thy(ii),1);
p = polyfit(ii*N_T,r2R(ii),1);
fprintf("px = %f, py = %f, p = %f\npthx = %f, pthy = %f",px(1), py(1), p(1), pthx(1), pthy(1));
%%   画出扩散系数和R的关系
D1 = [0.006446 0.064351 0.227403 1.550708];
D1 = D1/1e-4;
D2 = [0.000736 0.007454 0.018093 0.144077];
D2 = D2/1e-5;
R = [1.1 1.5 2 5];
x = (2^0.5.*R - 1)./(R - 1);%   圆阵列的实际孔喉比
figure;
hold on;
box on;
set(gca,'xscale','log')
yyaxis left;
plot(D1,R,'o-b','displayname','r=1e-4m');
plot(D2,R,'o-r','displayname','r=1e-5m');
xlabel('扩散系数/半径','fontsize',14);
ylabel('R——(圆直径+喉道宽度)/圆直径','fontsize',14);
title('修正扩散系数-R','fontsize',14);
yyaxis right;
plot(D1,x,'o--b','displayname','r=1e-4m');
plot(D2,x,'o--r','displayname','r=1e-5m');
ylabel('孔喉比','FontSize',14);
ax = gca;
yyaxis right;
ax.YDir = 'reverse';


function s = Sf(X1, X2, X3)
%   根据三角形三个顶点的坐标计算面积
a = norm(X2 - X1);
b = norm(X3 - X2);
c = norm(X2 - X3);
p = (a + b + c) / 2;
s = (p * (p - a) * (p - b) * (p - c))^0.5;
end





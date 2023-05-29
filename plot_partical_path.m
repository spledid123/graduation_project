%   画出单个粒子的某段轨迹以及周围的多孔介质
%%  画出圆阵列
clear;
figure;
axis off;
axis equal;
hold on;

r = 1;% 半径
R = [1.1 1.1];%    xy方向的圆心距/直径
%  xy方向的序号范围，圆心位置为r*Rx*(2*kx+1),r*Ry*(2*ky+1)
kx = 0:3;
ky = -10:1;
Rx = R(1);%  x/y方向圆心距/直径
Ry = R(2);

numx = length(kx);
numy = length(ky);
for i = 1:numx
    x0(i) = r * Rx * (2 * kx(i) + 1);
end
for i = 1:numy
    y0(i) = r * Ry * (2 * ky(i) + 1);
end
for i = 1:numx
    for j = 1:numy
        plot_cir(x0(i),y0(j),r,1);
    end
end

%%   画孔喉六边形模型
%   这样画没办法填充，只能借助ps
figure;
hold on;
axis equal;
axis off;

r = 5;% 半径
a = 0.5;%  喉道长度/直径
R = 1.2;%   （喉道宽度+直径)/直径
%   序号范围
kx = -1:1;
ky = -4:3;

t = (1 + a) * r;%   六边形边长的一半
l = 1.7320508075688772935274463415059;% 3^0.5
y = l * t * ky;
kx1 = 6 * kx;
kx2 = 6 * kx + 4;
kx3 = 6 * kx + 1;
kx4 = 6 * kx + 3;
for j = 1:length(ky)
    for i = 1:length(kx)
        if(mod(ky(j),2) == 0)
            plot_onecir_in_six(t * kx1(i), y(j), r, 2 * (R - 1) * r,2 * a * r, 0);
            plot_onecir_in_six(t * kx2(i), y(j), r, 2 * (R - 1) * r,2 * a * r, 1);
        else
            plot_onecir_in_six(t * kx3(i), y(j), r, 2 * (R - 1) * r,2 * a * r, 1);
            plot_onecir_in_six(t * kx4(i), y(j), r, 2 * (R - 1) * r,2 * a * r, 0);
        end
    end
end
%%  画出周期性圆阵列，九宫格示意图
clear;
figure;
axis off;
axis equal;
hold on;

%   x = 0左右两边的区块,圆阵列,如果一边是空白，numl = 0,如果区块的左下角(右下角)不是(0,0),调整tsdr/tsdl
numr = 50;    % 右边圆阵列圆数量
numl = 50;    % 左边圆阵列圆数量，没有记为0
Rr = 2;       % 右边圆阵列圆半径
Rl = 1.5;      % 左边圆阵列圆半径
tr = 2.2;      % 右边圆阵列圆半径+喉道宽度一半
tl = 2.2;      % 左边圆阵列圆半径+喉道宽度一半
tsdr = 0;      %  圆阵列所在正方形与区块的间隔
tsdl = 0;
if(numl ~= 0 && numr ~= 0 && numr * tr + tsdr ~= numl * tl + tsdl || tr < Rr || tl < Rl)
    error('参数设置错误');
end
%   右边区块
kx = 0:numr - 1;
ky = 0:numr - 1;
for i = 1:numr
    x0(i) = tr*(2*kx(i)+1) + tsdr;
end
for i = 1:numr
    y0(i) = tr*(2*ky(i)+1) + tsdr;
end
for i = 1:numr
    for j = 1:numr
        plot_cir(x0(i),y0(j),Rr,1);
    end
end
%   左边区块
kx = 0:numl - 1;
ky = 0:numl - 1;
for i = 1:numl
    x0(i) = -(tl*(2*kx(i)+1)) - tsdl;
end
for i = 1:numl
    y0(i) = tl*(2*ky(i)+1) + tsdl;
end
for i = 1:numl
    for j = 1:numl
        plot_cir(x0(i),y0(j),Rl,1);
    end
end
rx = 2 * (tr * numr + tsdr);
plot([-rx rx rx -rx -rx], [0 0 rx rx 0],'b--','LineWidth',2);
plot([0 0], [0 rx],'b--','LineWidth',2);
%%  画开放空间内圆阵列多孔介质模型
r = 200;    %   大圆半径
R = [2 2];%    xy方向的圆心距离/直径
xmax = 30;  %   每行圆个数
ymax = 30;  %   每列圆个数
kx = 0 : xmax - 1;
ky = 0 : ymax - 1;
rx = (xmax) * R(1);
ry = (ymax) * R(2);
for i = 1:xmax
    x0(i) = R(1) * (2 * kx(i) + 1) - rx;
end
for i = 1:ymax
    y0(i) = R(2) * (2 * ky(i) + 1) - ry;
end
figure;
hold on;
for i = 1:xmax
    for j = 1:ymax
        plot_cir(x0(i),y0(j),1,1);
    end
end
axis off;
axis equal;
plot_cir(0,0,r,'var','b--');
plot([-rx rx rx -rx -rx], [-ry -ry ry ry -ry],'b--','LineWidth',2);
%% 气体内随机圆多孔介质
clear;
r = 100;    %   大圆半径
rx = 25 * 2;    %   多孔介质范围
ry = 25 * 2;
R = [];     %   小圆半径
num = 150;  %   小圆个数
rmax = 5;   %   小圆半径范围
rmin = 0.1;
rx1 = rx - rmax;  %   圆心范围
ry1 = ry - rmax;
for i = 1:num
    x0(i) = rand() * 2 * rx1 - rx1;
    y0(i) = rand() * 2 * ry1 - ry1;
    R(i) = rand() * (rmax - rmin) + rmin;
end
hold on;
for i=1:num
    plot_cir(x0(i),y0(i),R(i),1);
end
axis off
axis equal
plot_cir(0,0,r,'var','b--');
plot([-rx rx rx -rx -rx], [-ry -ry ry ry -ry],'b--','LineWidth',2);
%%  气体内多孔介质模型——三角
s = [84 * 3^0.5, 1];
figure;
hold on;
axis equal;
axis off;
plot_rec([0, -105], s, 0, 1);
plot_rec([-52.5 * 3^0.5, 52.5], s, pi/3, 1);
plot_rec([52.5 * 3^0.5, 52.5], s, -pi/3, 1);
plot_cir(0,0,500,'var','b--');
%%  气体内多孔介质模型——四圆
figure;
hold on;
axis equal;
axis off;
plot_cir(0,0,500,'var','b--');
plot_cir(110,110,100,1);
plot_cir(110,-110,100,1);
plot_cir(-110,110,100,1);
plot_cir(-110,-110,100,1);
%%  气体内多孔介质模型——二圆
figure;
hold on;
axis equal;
axis off;
plot_cir(0,0,500,'var','b--');
plot_arc([-150, 0], 200, 205, 50/180*pi, 2*pi-50/180*pi);
plot_arc([150, 0], 200, 205, pi+50/180*pi, 3*pi-50/180*pi);
%%  气体内多孔介质模型——四正方形
figure;
hold on;
axis equal;
axis off;
plot_cir(0,0,500,'var','b--');
s = [50 50];
plot_rec([100 100],s,0,1);
plot_rec([100 -100],s,0,1);
plot_rec([-100 100],s,0,1);
plot_rec([-100 -100],s,0,1);
rx = 150;
plot([-rx -rx rx rx -rx], [rx -rx -rx rx rx],'b--','LineWidth',2);
%%  气体内多孔介质模型——正方形阵列
figure;
hold on;
axis equal;
axis off;
plot_cir(0,0,200,'var','b--');
s = [3 3];
xmax = 20;
ymax = 20;
kx = 0:xmax - 1;
ky = 0:ymax - 1;
R = s.*[1.5 1.5];%    xy方向的圆心距离/直径
rx = xmax * R(1);
ry = ymax * R(2);
for i = 1:xmax
    x0(i) = R(1) * (2 * kx(i) + 1) - rx;
end
for i = 1:ymax
    y0(i) = R(2) * (2 * ky(i) + 1) - ry;
end
for i = 1:xmax
    for j = 1:ymax
        plot_rec([x0(i),y0(j)],s,0,1);
    end
end
rx = R(1) * xmax - 3;
ry = R(2) * ymax - 3;
plot([-rx -rx rx rx -rx], [ry -ry -ry ry ry],'b--');
%%  气体内多孔介质模型——凹圆弧
figure;
hold on;
axis equal;
axis off;
r = [120 125];
plot_cir(0,0,500,'var','b--');
plot_arc([0, -150], r(1), r(2), pi/2+50/180*pi, pi/2-50/180*pi);
plot_arc([150/2*3^0.5, 150/2], r(1), r(2), 2*pi/3+pi/2+50/180*pi, 2*pi/3+pi/2-50/180*pi);
plot_arc([-150/2*3^0.5, 150/2], r(1), r(2), -2*pi/3+pi/2+50/180*pi, -2*pi/3+pi/2-50/180*pi);
%   画多孔介质边界
r0 = [0 -150];
r1 = [150/2*3^0.5, 150/2];
x1 = [r0(1) + 125*cos(40/180*pi) r0(2) + 125*sin(40/180*pi)];
x2 = [r1(1) + 125*cos(2*pi/3+pi/2+50/180*pi) r1(2) + 125*sin(2*pi/3+pi/2+50/180*pi)];
x3 = [r1(1) + 125*cos(2*pi/3+pi/2-50/180*pi) r1(2) + 125*sin(2*pi/3+pi/2-50/180*pi)];
x11 = [-x1(1) x1(2)];
x21 = [-x2(1) x2(2)];
x31 = [-x3(1) x3(2)];
plot([x1(1) x2(1)], [x1(2) x2(2)],'--b');
plot([x11(1) x21(1)], [x11(2) x21(2)],'--b');
plot([x3(1) x31(1)], [x3(2) x31(2)],'--b');
%   找"三角形"的三个顶点
X1 = [0 x2(2) + (-x2(1)) * (x2(2) - x1(2)) / (x2(1) - x1(1))];
X2 = [(x3(2) - x2(2)) * (x2(1) - x1(1)) / (x2(2) - x1(2)) + x2(1) x3(2)];
X3 = [-X2(1) X2(2)];
plot([X1(1) X2(1) X3(1) X1(1)], [X1(2) X2(2) X3(2) X1(2)],'--b','LineWidth',2);
%%  气体内多孔介质模型——凸形弧
figure;
hold on;
axis equal;
axis off;
r = [150 155];
plot_cir(0,0,500,'var','b--');
plot_arc([0, -200], r(1), r(2), pi/2-30/180*pi, pi/2+30/180*pi);
plot_arc([200/2*3^0.5, 200/2], r(1), r(2), 2*pi/3+pi/2-30/180*pi, 2*pi/3+pi/2+30/180*pi);
plot_arc([-200/2*3^0.5, 200/2], r(1), r(2), -2*pi/3+pi/2-30/180*pi, -2*pi/3+pi/2+30/180*pi);
%%  气体内多孔介质模型——开口圆弧
figure;
hold on;
axis equal;
axis off;
r = [100 100+1];
plot_cir(0,0,200,'var','b--');
plot_cir(0,0,100,'var','b--');
plot_arc([0, 0], r(1), r(2), -pi/6+5/180*pi, pi/2-5/180*pi);
plot_arc([0, 0], r(1), r(2), pi/2+5/180*pi, -5*pi/6-5/180*pi);
plot_arc([0, 0], r(1), r(2), -5*pi/6+5/180*pi, -pi/6-5/180*pi);
%%  随机圆阵列
figure;
hold on;
axis equal;
axis off;
cir = load('cir.mat');
cir = cir.cir_5;
len = size(cir);
len = len(1) - 1;
for i = 1:len
    plot_cir(cir(i,1), cir(i,2), cir(i,3), 1);
end
r = 500;
plot_cir(0, 0, 500, 'var', 'b--');
rx = 205;
% plot([-rx -rx rx rx -rx], [rx -rx -rx rx rx],'b--','LineWidth',2);
xlim([-205 205]);
ylim([-205 205]);
%%  随机矩形
figure;
hold on;
axis equal;
axis off;
cir = load('cube.mat');
cir = cir.cube_1;
len = size(cir);
len = len(1) - 1;
for i = 1:len
    plot_rec([cir(i,1), cir(i,2)], [cir(i,3), cir(i,4)], cir(i,5), 1);
end
r = 500;
plot_cir(0, 0, 500, 'var', 'b--');
rx = 205;
% plot([-rx -rx rx rx -rx], [rx -rx -rx rx rx],'b--','LineWidth',2);
xlim([-205 205]);
ylim([-205 205]);
%%  画单粒子轨迹
% clc;
clear;
filename = "data\path_cubecir\path_2000_R_1.1\rxT_circle_onepartical_1_.txt";
A = readtable(filename);
Tab = 0;%   平均吸附时间
T = get_T(filename,Tab);
len = size(A);
len = len(1);
%   轨迹范围
path_start = 1;
path_end = len;
rx=table2array(A(path_start:path_end,2));
ry=table2array(A(path_start:path_end,3));
T = T(path_start:path_end);
ry(end)=NaN;
T(end)=NaN;
%   不带颜色
hold on;
axis equal;
plot(rx,ry);
box on;
%   根据时间附上颜色
% figure;
% hold on;
% axis equal;
% patch(rx,ry,T,'EdgeColor','interp','Marker','.','MarkerFaceColor','flat');
% c1=colorbar;
% set(get(c1,'title'),'string','时间','fontsize',14);
% box on;



function T = get_T(filename,T0)
%   根据单粒子轨迹文件，返回每个点的对应时间
%   T0为平均吸附时间
A = readtable(filename);
X = table2array(A(:,2));%  x向量
Y = table2array(A(:,3));%  y向量
%   每一步步长
DX = [];    %   每一步走的X分量
DY = [];    %   每一步走的Y分量
DR = [];     %   每一步走的长度

len = length(X);
for i = 2:len
    DX(i-1) = X(i) - X(i - 1);
    DY(i-1) = Y(i) - Y(i - 1);
    DR=(DX(i-1)^2+DY(i-1)^2)^0.5;
end
DX = DX(1:len-1);
DY = DY(1:len-1);
DR = DX.^2 + DY.^2;
DR = DR.^0.5;
%   选择轨迹段
T = [];
T(1) = 0;
filename_v = 'data\v_10000\Blotzmann_water_T_40_2D.txt';
Av = readtable(filename_v);
Av = table2array(Av);
for i = 2:len
    rv = randi(10000);
    rv = Av(rv);
    Tab = get_T_ab(T0);
    T(i) = T(i - 1) + DR(i - 1) / 1e5 / rv + Tab;
end
end
function s = get_T_ab(T)
%   给定平均时间，得到一个指数分布的吸附时间
    if(T == 0)
        s = 0;
        return;
    end
    s = -T*log(rand());
end
function []=plot_onecir_in_six(x,y,r,h,len,flag)
%   给出圆心和类型,画出对应孔喉
%   [x,y]为圆心,r为半径,h为孔的半径,len为孔的长度
%   flag==0表示为(0,0)处的孔;flag==1表示为(-2t,0)处的孔

if(flag==0)
    k = (r^2-(r-h/2)^2)/2/r/(r-h/2);
    th = pi+atan(k):0.01:4*pi/3;
    plot_pore_cir(x,y,r,th);
    th = th - 2 * pi / 3;%  顺时针旋转120°
    plot_pore_cir(x,y,r,th);
    th = th -2 * pi / 3;
    plot_pore_cir(x,y,r,th);
    th = 2*pi/3:0.001:pi-atan(k);
    plot_pore_cir(x,y,r,th);
    th = th - 2*pi/3;
    plot_pore_cir(x,y,r,th);
    th = th - 2*pi/3;
    plot_pore_cir(x,y,r,th);

    rx = [-r;-k*r];
    rr = -h/2-rx(2);
    th = pi/2:-0.001:atan(k);
    plot_pore_cir(x+rx(1),y+rx(2),rr,th);
    rx = [-0.5 3^0.5/2;-3^0.5/2 -0.5]*rx;
    th = th - 2*pi/3;
    plot_pore_cir(x+rx(1),y+rx(2),rr,th);
    th = th - 2*pi/3;
    rx = [-0.5 3^0.5/2;-3^0.5/2 -0.5]*rx;
    plot_pore_cir(x+rx(1),y+rx(2),rr,th);
    rx = [-r;k*r];
    th = -pi/2:0.001:-atan(k);
    plot_pore_cir(x+rx(1),y+rx(2),rr,th);
    rx = [-0.5 3^0.5/2;-3^0.5/2 -0.5]*rx;
    th = th - 2*pi/3;
    plot_pore_cir(x+rx(1),y+rx(2),rr,th);
    th = th - 2*pi/3;
    rx = [-0.5 3^0.5/2;-3^0.5/2 -0.5]*rx;
    plot_pore_cir(x+rx(1),y+rx(2),rr,th);
 
    rx1 = [-r;h/2];
    rx2 = [-r;-h/2];
    rx3 = [-r - len/2; h/2];
    rx4 = [-r - len/2; -h/2];
    line([x + rx1(1) x + rx3(1)],[y + rx1(2) y + rx3(2)]);
    line([x + rx2(1) x + rx4(1)],[y + rx2(2) y + rx4(2)]);
    rx1 = [-0.5 3^0.5/2;-3^0.5/2 -0.5]*rx1;
    rx2 = [-0.5 3^0.5/2;-3^0.5/2 -0.5]*rx2;
    rx3 = [-0.5 3^0.5/2;-3^0.5/2 -0.5]*rx3;
    rx4 = [-0.5 3^0.5/2;-3^0.5/2 -0.5]*rx4;
    line([x + rx1(1) x + rx3(1)],[y + rx1(2) y + rx3(2)]);
    line([x + rx2(1) x + rx4(1)],[y + rx2(2) y + rx4(2)]);
    rx1 = [-0.5 3^0.5/2;-3^0.5/2 -0.5]*rx1;
    rx2 = [-0.5 3^0.5/2;-3^0.5/2 -0.5]*rx2;
    rx3 = [-0.5 3^0.5/2;-3^0.5/2 -0.5]*rx3;
    rx4 = [-0.5 3^0.5/2;-3^0.5/2 -0.5]*rx4;
    line([x + rx1(1) x + rx3(1)],[y + rx1(2) y + rx3(2)]);
    line([x + rx2(1) x + rx4(1)],[y + rx2(2) y + rx4(2)]);
else
    k = (r^2-(r-h/2)^2)/2/r/(r-h/2);
    th = -atan(k):-0.01:-pi/3;
    plot_pore_cir(x,y,r,th);
    th = th - 2 * pi / 3;%  顺时针旋转120°
    plot_pore_cir(x,y,r,th);
    th = th -2 * pi / 3;
    plot_pore_cir(x,y,r,th);
    th = pi/3:-0.001:atan(k);
    plot_pore_cir(x,y,r,th);
    th = th - 2*pi/3;
    plot_pore_cir(x,y,r,th);
    th = th - 2*pi/3;
    plot_pore_cir(x,y,r,th);

    rx = [r;-k*r];
    rr = -h/2-rx(2);
    th = pi/2:0.001:pi-atan(k);
    plot_pore_cir(x+rx(1),y+rx(2),rr,th);
    rx = [-0.5 3^0.5/2;-3^0.5/2 -0.5]*rx;
    th = th - 2*pi/3;
    plot_pore_cir(x+rx(1),y+rx(2),rr,th);
    th = th - 2*pi/3;
    rx = [-0.5 3^0.5/2;-3^0.5/2 -0.5]*rx;
    plot_pore_cir(x+rx(1),y+rx(2),rr,th);
    rx = [r;k*r];
    th = -pi/2:-0.001:-pi+atan(k);
    plot_pore_cir(x+rx(1),y+rx(2),rr,th);
    rx = [-0.5 3^0.5/2;-3^0.5/2 -0.5]*rx;
    th = th - 2*pi/3;
    plot_pore_cir(x+rx(1),y+rx(2),rr,th);
    th = th - 2*pi/3;
    rx = [-0.5 3^0.5/2;-3^0.5/2 -0.5]*rx;
    plot_pore_cir(x+rx(1),y+rx(2),rr,th);
 
    rx1 = [r;h/2];
    rx2 = [r;-h/2];
    rx3 = [r + len/2; h/2];
    rx4 = [r + len/2; -h/2];
    line([x + rx1(1) x + rx3(1)],[y + rx1(2) y + rx3(2)]);
    line([x + rx2(1) x + rx4(1)],[y + rx2(2) y + rx4(2)]);
    rx1 = [-0.5 3^0.5/2;-3^0.5/2 -0.5]*rx1;
    rx2 = [-0.5 3^0.5/2;-3^0.5/2 -0.5]*rx2;
    rx3 = [-0.5 3^0.5/2;-3^0.5/2 -0.5]*rx3;
    rx4 = [-0.5 3^0.5/2;-3^0.5/2 -0.5]*rx4;
    line([x + rx1(1) x + rx3(1)],[y + rx1(2) y + rx3(2)]);
    line([x + rx2(1) x + rx4(1)],[y + rx2(2) y + rx4(2)]);
    rx1 = [-0.5 3^0.5/2;-3^0.5/2 -0.5]*rx1;
    rx2 = [-0.5 3^0.5/2;-3^0.5/2 -0.5]*rx2;
    rx3 = [-0.5 3^0.5/2;-3^0.5/2 -0.5]*rx3;
    rx4 = [-0.5 3^0.5/2;-3^0.5/2 -0.5]*rx4;
    line([x + rx1(1) x + rx3(1)],[y + rx1(2) y + rx3(2)]);
    line([x + rx2(1) x + rx4(1)],[y + rx2(2) y + rx4(2)]);
end
end
function []=plot_pore_cir(x,y,r,th)
%   画圆和临近孔，只画一部分，th给出圆的范围
    rx = x + r * cos(th);
    ry = y + r * sin(th);
    plot(rx,ry,'r');
end
function [] = plot_arc(c, r1, r2, th1, th2,varargin)
%   画圆弧（圆环的部分）
%   c=[cx,cy]为圆心，r1<r2为圆环半径，th1-th2逆时针为圆弧角度
%   可选参数flag默认为1，表示填充，flag == 0不填充
%   名称-数值对:var为plot的参数；fillvar为填充颜色
p = inputParser;            % 函数的输入解析器
fcn = @(x) x == 0 || x == 1;
addOptional(p,'flag',1,fcn);
addParameter(p,'var',"r-");   
addParameter(p,'fillvar',"r");      % 默认参数
parse(p,varargin{:});
th1 = th1 - 2*pi*floor(th1/2/pi);
th2 = th2 - 2*pi*floor(th2/2/pi);
if(th1 > th2)
    th1 = th1 - 2*pi;
end
th = [th1:0.01:th2 th2];
x1 = c(1) + r1 * cos(th);
y1 = c(2) + r1 * sin(th);
x2 = c(1) + r2 * cos(flip(th));
y2 = c(2) + r2 * sin(flip(th));
plot([x1 x2 x1(1)], [y1 y2 y1(1)],p.Results.var);
if(p.Results.flag)
    fill([x1 x2 x1(1)], [y1 y2 y1(1)],p.Results.fillvar);
end
end
%%  三外圆切圆
clc;
clear;
th = 5;%    开口角度的一半
th = 30 + th;
b = (4 - 2 * sin(th/180*pi));
c = 2 - 4 * sin(th/180*pi);
a = (-b + (b^2 - 4*c)^0.5)/2 * 200;
%   画出四个圆,一个圆被其它三个切，留下开口为10°
figure;plot_cir(0,0,200,1);hold on;axis equal;axis off;xlim([-205 205]);ylim([-205 205]);
set(gca,'position',[0,0,1,1]);
set(gcf,'position',[0,0,800,800]);
saveas(gcf,'data\pic\cir1.bmp');
close();
figure;plot_cir(0,-2*200-a,3^0.5*200,1);hold on;axis equal;axis off;xlim([-205 205]);ylim([-205 205]);
set(gca,'position',[0,0,1,1]);
set(gcf,'position',[0,0,800,800]);
saveas(gcf,'data\pic\cir2.bmp');
close();
figure;plot_cir(-3^0.5*200-cos(30/180*pi)*a,200+a/2,3^0.5*200,1);hold on;axis equal;axis off;xlim([-205 205]);ylim([-205 205]);
set(gca,'position',[0,0,1,1]);
set(gcf,'position',[0,0,800,800]);
saveas(gcf,'data\pic\cir3.bmp');
close();
figure;plot_cir(3^0.5*200+cos(30/180*pi)*a,200+a/2,3^0.5*200,1);hold on;axis equal;axis off;xlim([-205 205]);ylim([-205 205]);
set(gca,'position',[0,0,1,1]);
set(gcf,'position',[0,0,800,800]);
saveas(gcf,'data\pic\cir4.bmp');
close();
%   图像计算操作
cir1 = imread('data\pic\cir1.bmp');
cir2 = imread('data\pic\cir2.bmp');
cir3 = imread('data\pic\cir3.bmp');
cir4 = imread('data\pic\cir4.bmp');
cir_1 = red_lap(cir1, cir2);
cir_2 = red_lap(cir1, cir3);
cir_3 = red_lap(cir1, cir4);
cir = red_add(red_add(cir_1, cir_2), cir_3);
figure;imshow(cir);

%%   圆内三角
clc;
clear;
th = 5;%    开口圆心角的一半
a = 200 * (sin(th/180*pi) * 3^0.5 + cos(th/180*pi));
figure;plot_cir(0,0,200,1);hold on;axis equal;axis off;xlim([-205 205]);ylim([-205 205]);
set(gca,'position',[0,0,1,1]);
set(gcf,'position',[0,0,800,800]);
saveas(gcf,'data\pic\cir1.bmp');
close();
figure;fill([0 -cos(30/180*pi)*a cos(30/180*pi)*a], [-a sin(30/180*pi)*a sin(30/180*pi)*a],'r','LineStyle','none');
hold on;axis equal;axis off;xlim([-205 205]);ylim([-205 205]);
set(gca,'position',[0,0,1,1]);
set(gcf,'position',[0,0,800,800]);
saveas(gcf,'data\pic\tri.bmp');
close();
cir1 = imread('data\pic\cir1.bmp');
tri = imread('data\pic\tri.bmp');
cir = red_minus(cir1, tri);
figure;
imshow(cir);

%%  三角内直边
clc;
clear;
t = 1/10;%   开口占边长的比例
th = 0/180*pi;%    逆时针旋转角
plot_tri_rec(t, 200, [100 0], th, 205, 'data\pic\tr.bmp');
%%  三角内凸边
clc;
clear;
t = 1/10;%   开口占边长的比例
th = 0;%    逆时针旋转角
plot_tri_vex(t, 200, [0 0], th, 205, 'data\pic\tr.bmp');
%%  三角内凹边
clc;
clear;
t = 1/50;%   开口占边长的比例
th = 0;%    逆时针旋转角
plot_tri_cave(t, 200, [0 0], th, 205, 'data\pic\tr.bmp');
%%  三角形正六边形阵列
clc;
clear;
t = 1/10;%   开口占边长的比例
a = 200;%   边长
%   六个三角形的信息
c1 = [a/2 a/2/3^0.5];
c2 = [0 a/3^0.5];
c3 = [-a/2 a/2/3^0.5];
c4 = [-a/2 -a/2/3^0.5];
c5 = [0 -a/3^0.5];
c6 = [a/2 -a/2/3^0.5];
plot_tri_cave(t, a, c1, 0, 205, 'data\pic\1.bmp');close();
plot_tri_cave(t, a, c2, pi/3, 205, 'data\pic\2.bmp');close();
plot_tri_cave(t, a, c3, 0, 205, 'data\pic\3.bmp');close();
plot_tri_cave(t, a, c4, pi/3, 205, 'data\pic\4.bmp');close();
plot_tri_cave(t, a, c5, 0, 205, 'data\pic\5.bmp');close();
plot_tri_cave(t, a, c6, pi/3, 205, 'data\pic\6.bmp');close();
cir1 = imread('data\pic\1.bmp');
cir2 = imread('data\pic\2.bmp');
cir3 = imread('data\pic\3.bmp');
cir4 = imread('data\pic\4.bmp');
cir5 = imread('data\pic\5.bmp');
cir6 = imread('data\pic\6.bmp');
cir = red_add(cir1, red_add(cir2, red_add(cir3, red_add(cir4, red_add(cir5, cir6)))));
imshow(cir);
%%  双重三角形六边形阵列
clc;
clear;
t = 1/10;%   开口占边长的比例
a = 100;%   边长
tx = a/2;
ty = a/3^0.5/2;
%   六个三角形的信息
c1 = [a/2 a/2/3^0.5];
c2 = [0 a/3^0.5];
c3 = [-a/2 a/2/3^0.5];
c4 = [-a/2 -a/2/3^0.5];
c5 = [0 -a/3^0.5];
c6 = [a/2 -a/2/3^0.5];
c1x = [tx -tx 0 2*tx 0 -2*tx -3*tx 3*tx -2*tx 2*tx -tx tx];
c1y = [ty ty -2*ty 4*ty 4*ty 4*ty ty ty -2*ty -2*ty -5*ty -5*ty];
c2y = -c1y;
for i = 1:12
    plot_tri_rec(t, a, [c1x(i) c1y(i)], 0, 205, strcat('data\pic\', int2str(i), '.bmp'));
    close();
    plot_tri_rec(t, a, [c1x(i) c2y(i)], pi/3, 205, strcat('data\pic\', int2str(i + 12), '.bmp'));
    close();
end
cir = imread('data\pic\1.bmp');
for i = 2:24
    cir = red_add(cir, imread(strcat('data\pic\', int2str(i),'.bmp')));
end
imshow(cir);
%%  得到SDF
[s, sdf, xlim, ylim] = get_dist(gca,gcf);
close();
ppm = struct('name', 'tri', 'exp', '三角形凹圆,边长200,t=1/10,中心[100 0]', 'maxx', xlim, 'maxy', ylim, 'bina', s, 'sdf', sdf);
a = dir('data/porousmedia/*.mat');
len = length(a);len=13;
eval(strcat('pm', int2str(len + 1), '{1} = ppm;'));
figure;
eval(strcat('imshow(pm', int2str(len + 1), '{1}.bina);'));
%%  存储
a = dir('data/porousmedia/*.mat');
len = length(a);len=13;
filename = strcat('data\porousmedia\pm', int2str(len + 1), '.mat');
varname = strcat('pm', int2str(len + 1));
save(filename, varname,'-mat');
set(gcf,'position',[0 0 800 800]);
set(gca,'position',[0 0 1 1]);
saveas(gca,strcat('data\porousmedia\pm', int2str(len + 1), '.png'));
close();
%%  算孔隙度,Sp
t = 3;
lim = 205;
S = 100*300*3^0.5;
load(strcat('data\porousmedia\pm', int2str(t), '.mat'));
eval(strcat('ppm = pm', int2str(t), ';'));
phi = sum(ppm{1}.bina,'all')/numel(ppm{1}.bina);
Ss = lim^2 * (1 - phi);
Sp = S - Ss;
vpa(Sp);
%%
t = 14;
load(strcat('data\porousmedia\pm', int2str(t), '.mat'));
figure;
eval(strcat('imshow(pm', num2str(t), '{1}.bina)'));
set(gca, 'position', [0 0 1 1]);
set(gcf, 'position', [0 0 800 800]);
saveas(gca, strcat('data\porousmedia\pm', int2str(t), '.png'));
close();


function cir = red_lap(cir1, cir2)
%   用红色表现多孔介质，白色表示空间
%   cir1,cir2分别为rgb矩阵表示图形
%   cir表示cir1和cir2重叠的部分
cir1 = 255 - cir1;
cir2 = 255 - cir2;
cir = imsubtract(cir1,cir2);
cir1 = 255 - cir1;
cir = 255 - cir;
cir = imsubtract(cir,cir1);
cir = 255 - cir;
end
function cir = red_add(cir1, cir2)
%   用红色表现多孔介质，白色表示空间
%   cir1,cir2分别为rgb矩阵表示图形
%   cir表示cir1和cir2红色加起来
cir1 = 255 - cir1;
cir2 = 255 - cir2;
cir = imadd(cir1, cir2);
cir = 255 - cir;
end
function cir = red_minus(cir1, cir2)
%   用红色表现多孔介质，白色表示空间
%   cir1,cir2分别为rgb矩阵表示图形
%   cir表示cir1红色部分减去cir1,cir2红色重叠部分
cir2 = red_rev(cir2);
cir2 = 255 - cir2;
cir1 = 255 - cir1;
cir = imsubtract(cir1,cir2);
cir = 255 - cir;
cir1 = 255 - cir1;
cir = imsubtract(cir,cir1);
cir = 255 - cir;
end
function cir = red_rev(cir1)
%   用红色表现多孔介质，白色表示空间
%   cir1为rgb矩阵表示图形
%   cir表示cir1红色部分与白色部分调换
cir(:,:,2:3) = 255 -  cir1(:,:,2:3);
end
function [] = plot_tri_rec(t, a, c, th, lim, filename)
%   生成三角形内直边的图
%   t为开口比例,a为边长,c为中心点位置,th为逆时针旋转角,th = 0一条中线y轴
%   图像范围[-lim lim],[-lim lim]
figure;
r1 = [0 a/3^0.5] * [cos(th), sin(th); -sin(th) cos(th)];
r1 = r1 + c;
r2 = [-a/2 -a/3^0.5/2] * [cos(th), sin(th); -sin(th) cos(th)];
r2 = r2 + c;
r3 = [a/2 -a/3^0.5/2] * [cos(th), sin(th); -sin(th) cos(th)];
r3 = r3 + c;
plot_tri(r1, r2, r3, 1);
hold on;axis equal;axis off;xlim([-lim lim]);ylim([-lim lim]);
set(gca,'position',[0,0,1,1]);
set(gcf,'position',[0,0,800,800]);
saveas(gcf,'data\pic\tri.bmp');
close();
sd = lim/2050;% 矩形偏置以消除边缘
r4 = [-(a/4/3^0.5 + sd)/2*3^0.5 (a/4/3^0.5 + sd)/2] * [cos(th), sin(th); -sin(th) cos(th)];
r4 = c + r4;
figure;plot_rec(r4, [a/4/3^0.5 + sd a/2*t], -pi/6 + th, 1);
hold on;axis equal;axis off;xlim([-lim lim]);ylim([-lim lim]);
set(gca,'position',[0,0,1,1]);
set(gcf,'position',[0,0,800,800]);
saveas(gcf,'data\pic\rec1.bmp');
close();
r4 = [(a/4/3^0.5 + sd)/2*3^0.5 (a/4/3^0.5 + sd)/2] * [cos(th), sin(th); -sin(th) cos(th)];
r4 = c + r4;
figure;plot_rec(r4, [a/4/3^0.5 + sd a/2*t], pi/6 + th, 1);
hold on;axis equal;axis off;xlim([-lim lim]);ylim([-lim lim]);
set(gca,'position',[0,0,1,1]);
set(gcf,'position',[0,0,800,800]);
saveas(gcf,'data\pic\rec2.bmp');
close();
r4 = [0 -(a/4/3^0.5 + sd)] * [cos(th), sin(th); -sin(th) cos(th)];
r4 = c + r4;
figure;plot_rec(r4, [a/4/3^0.5 + sd a/2*t], pi/2 + th, 1);
hold on;axis equal;axis off;xlim([-lim lim]);ylim([-lim lim]);
set(gca,'position',[0,0,1,1]);
set(gcf,'position',[0,0,800,800]);
saveas(gcf,'data\pic\rec3.bmp');
close();
cir1 = imread('data\pic\tri.bmp');
cir2 = imread('data\pic\rec1.bmp');
cir3 = imread('data\pic\rec2.bmp');
cir4 = imread('data\pic\rec3.bmp');
cir = red_minus(cir1, cir2);
cir = red_minus(cir, cir3);
cir = red_minus(cir, cir4);
figure;
imshow(cir);
set(gca,'position',[0,0,1,1]);
set(gcf,'position',[0,0,800,800]);
saveas(gca, filename);
end

function [] = plot_tri_vex(t, a, c, th, lim, filename)
%   生成三角形内凸圆的图
%   t为开口比例,a为边长,c为中心点位置,th为逆时针旋转角,th = 0一条中线y轴
%   图像范围[-lim lim],[-lim lim]
figure;
r1 = [0 a/3^0.5] * [cos(th), sin(th); -sin(th) cos(th)];
r1 = r1 + c;
r2 = [-a/2 -a/3^0.5/2] * [cos(th), sin(th); -sin(th) cos(th)];
r2 = r2 + c;
r3 = [a/2 -a/3^0.5/2] * [cos(th), sin(th); -sin(th) cos(th)];
r3 = r3 + c;
plot_tri(r1, r2, r3, 1);
hold on;axis equal;axis off;xlim([-lim lim]);ylim([-lim lim]);
set(gca,'position',[0,0,1,1]);
set(gcf,'position',[0,0,800,800]);
saveas(gcf,'data\pic\tri.bmp');
close();
figure;plot_cir(r1(1), r1(2), a/2*(1-t), 1);
hold on;axis equal;axis off;xlim([-lim lim]);ylim([-lim lim]);
set(gca,'position',[0,0,1,1]);
set(gcf,'position',[0,0,800,800]);
saveas(gcf,'data\pic\cir1.bmp');
close();
figure;plot_cir(r2(1), r2(2), a/2*(1-t), 1);
hold on;axis equal;axis off;xlim([-lim lim]);ylim([-lim lim]);
set(gca,'position',[0,0,1,1]);
set(gcf,'position',[0,0,800,800]);
saveas(gcf,'data\pic\cir2.bmp');
close();
figure;plot_cir(r3(1), r3(2), a/2*(1-t), 1)
hold on;axis equal;axis off;xlim([-lim lim]);ylim([-lim lim]);
set(gca,'position',[0,0,1,1]);
set(gcf,'position',[0,0,800,800]);
saveas(gcf,'data\pic\cir3.bmp');
close();
cir1 = imread('data\pic\tri.bmp');
cir2 = imread('data\pic\cir1.bmp');
cir3 = imread('data\pic\cir2.bmp');
cir4 = imread('data\pic\cir3.bmp');
cir_1 = red_lap(cir1, cir2);
cir_2 = red_lap(cir1, cir3);
cir_3 = red_lap(cir1, cir4);
cir = red_add(red_add(cir_1, cir_2), cir_3);
figure;
imshow(cir);
set(gca,'position',[0,0,1,1]);
set(gcf,'position',[0,0,800,800]);
saveas(gcf, filename);
end

function [] = plot_tri_cave(t, a, c, th, lim, filename)
%   生成三角形内凹圆的图
%   t为开口比例,a为边长,c为中心点位置,th为逆时针旋转角,th = 0一条中线y轴
%   图像范围[-lim lim],[-lim lim]
figure;
r1 = [0 a/3^0.5] * [cos(th), sin(th); -sin(th) cos(th)];
r1 = r1 + c;
r2 = [-a/2 -a/3^0.5/2] * [cos(th), sin(th); -sin(th) cos(th)];
r2 = r2 + c;
r3 = [a/2 -a/3^0.5/2] * [cos(th), sin(th); -sin(th) cos(th)];
r3 = r3 + c;
plot_tri(r1, r2, r3, 1);
hold on;axis equal;axis off;xlim([-lim lim]);ylim([-lim lim]);
set(gca,'position',[0,0,1,1]);
set(gcf,'position',[0,0,800,800]);
saveas(gcf,'data\pic\tri.bmp');
close();
figure;
plot_cir(c(1), c(2), ((a/3^0.5/2)^2 + (a/3^0.5/2)^2*3*t^2)^0.5, 1);
hold on;axis equal;axis off;xlim([-lim lim]);ylim([-lim lim]);
set(gca,'position',[0,0,1,1]);
set(gcf,'position',[0,0,800,800]);
saveas(gcf,'data\pic\cir1.bmp');
close();
cir1 = imread('data\pic\tri.bmp');
cir2 = imread('data\pic\cir1.bmp');
cir = red_minus(cir1, cir2);
figure;
imshow(cir);
set(gca,'position',[0,0,1,1]);
set(gcf,'position',[0,0,800,800]);
saveas(gcf, filename);
end
function s = Shexagon(a)
%   正六边形的面积，a为边长
s = 3^1.5/2*a*a;
end

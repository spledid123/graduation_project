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
figure;fill([0 -cos(30/180*pi)*a cos(30/180*pi)*a], [-a sin(30/180)*pi*a sin(30/180)*pi*a],'r','LineStyle','none');
xlim([-205 205]);ylim([-205 205]);
set(gca,'position',[0,0,1,1]);
set(gcf,'position',[0,0,800,800]);
saveas(gcf,'data\pic\tri.bmp');
close();
cir1 = imread('data\pic\cir1.bmp');
tri = imread('data\pic\tri.bmp');
cir = red_minus(cir1, tri);
imshow(cir);

%%  三角内直边
clc;
clear;
t = 1/10;%   开口占边长的比例
figure;
fill([0 100*3^0.5 -100*3^0.5], [200 -100 -100],'r','LineStyle','none');
hold on;axis equal;axis off;xlim([-205 205]);ylim([-205 205]);
set(gca,'position',[0,0,1,1]);
set(gcf,'position',[0,0,800,800]);
saveas(gcf,'data\pic\tri.bmp');
close();
sd = 1;% 矩形偏置以消除边缘
figure;plot_rec([-25*3^0.5 - sd*3^0.5 25 + sd], [50 + 2*sd 100*t*3^0.5], -pi/6, 1);
hold on;axis equal;axis off;xlim([-205 205]);ylim([-205 205]);
set(gca,'position',[0,0,1,1]);
set(gcf,'position',[0,0,800,800]);
saveas(gcf,'data\pic\rec1.bmp');
close();
figure;plot_rec([25*3^0.5 + sd*3^0.5 25 + sd], [50 + 2*sd 100*t*3^0.5], pi/6, 1);
hold on;axis equal;axis off;xlim([-205 205]);ylim([-205 205]);
set(gca,'position',[0,0,1,1]);
set(gcf,'position',[0,0,800,800]);
saveas(gcf,'data\pic\rec2.bmp');
close();
figure;plot_rec([0 -50 - sd*2], [50 + 2*sd 100*t*3^0.5], pi/2, 1);
hold on;axis equal;axis off;xlim([-205 205]);ylim([-205 205]);
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
imshow(cir);
%%  三角内凸边
clc;
clear;
t = 1/10;%   开口占边长的比例
figure;
fill([0 100*3^0.5 -100*3^0.5], [200 -100 -100],'r','LineStyle','none');
hold on;axis equal;axis off;xlim([-205 205]);ylim([-205 205]);
set(gca,'position',[0,0,1,1]);
set(gcf,'position',[0,0,800,800]);
saveas(gcf,'data\pic\tri.bmp');
close();
figure;plot_cir(0, 200, 100*3^0.5*(1-t), 1);
hold on;axis equal;axis off;xlim([-205 205]);ylim([-205 205]);
set(gca,'position',[0,0,1,1]);
set(gcf,'position',[0,0,800,800]);
saveas(gcf,'data\pic\cir1.bmp');
close();
figure;plot_cir(100*3^0.5, -100, 100*3^0.5*(1-t), 1);
hold on;axis equal;axis off;xlim([-205 205]);ylim([-205 205]);
set(gca,'position',[0,0,1,1]);
set(gcf,'position',[0,0,800,800]);
saveas(gcf,'data\pic\cir2.bmp');
close();
figure;plot_cir(-100*3^0.5, -100, 100*3^0.5*(1-t), 1)
hold on;axis equal;axis off;xlim([-205 205]);ylim([-205 205]);
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
imshow(cir);
%%  三角内凹边
clc;
clear;
t = 1/10;%   开口占边长的比例
figure;
fill([0 100*3^0.5 -100*3^0.5], [200 -100 -100],'r','LineStyle','none');
hold on;axis equal;axis off;xlim([-205 205]);ylim([-205 205]);
set(gca,'position',[0,0,1,1]);
set(gcf,'position',[0,0,800,800]);
saveas(gcf,'data\pic\tri.bmp');
close();
figure;
plot_cir(0,0,(100^2 + 100^2*3*t^2)^0.5, 1);
hold on;axis equal;axis off;xlim([-205 205]);ylim([-205 205]);
set(gca,'position',[0,0,1,1]);
set(gcf,'position',[0,0,800,800]);
saveas(gcf,'data\pic\cir1.bmp');
close();
cir1 = imread('data\pic\tri.bmp');
cir2 = imread('data\pic\cir1.bmp');
cir = red_minus(cir1, cir2);
imshow(cir);

%%  得到SDF
[s, sdf, xlim, ylim] = get_dist(gca,gcf);
ppm = struct('name', 'tri3', 'exp', '单体为三角形，被凹圆切', 'maxx', xlim, 'maxy', ylim, 'bina', s, 'sdf', sdf);
load('data\pm.mat');
pm{length(pm) + 1} = ppm;
save('data\pm.mat', 'pm');


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


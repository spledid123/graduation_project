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
cir1 = 255 - cir1;
cir2 = 255 - cir2;
cir3 = 255 - cir3;
cir4 = 255 - cir4;
cir = imsubtract(cir1,cir2);
cir = imsubtract(cir,cir3);
cir = imsubtract(cir,cir4);
cir1 = 255 - cir1;
cir = 255 - cir;
cir = imsubtract(cir,cir1);
cir = 255 - cir;
imshow(cir);

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
tri(:,:,2:3) = 255 -  tri(:,:,2:3);
cir1 = 255 - cir1;
tri = 255 - tri;
cir = imsubtract(cir1,tri);
cir1 = 255 - cir1;
cir = 255 - cir;
cir = imsubtract(cir,cir1);
cir = 255 - cir;
imshow(cir);


%%  得到SDF
[s, sdf, xlim, ylim] = get_dist(gca,gcf);
pm{2} = struct('name', tri, 'exp', '单体为圆，被三角形切', 'maxx', xlim, 'maxy', ylim, 'bina', s, 'sdf', sdf);



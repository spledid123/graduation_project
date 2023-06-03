function [s, sdf, sumx, sumy] = get_dist(gca, gcf)
%   由rgb图像转为二值图像s，超分辨率，得到xy方向的像素个数sumx,sumy（固定为1e4),在算出SDF
%   SDF由bwdist得到，每个像素的大小为1，固体内为负，外为正
%   要求固体与空白为不同颜色，bulk为白色
set(gca,'position',[0,0,1,1]);
set(gcf,'position',[0,0,800,800]);
axis equal;
saveas(gca,'data\pic\get_dist_pic.bmp');
t = imread('data\pic\get_dist_pic.bmp');
ts = rgb2gray(t);%  灰度图
%   插值超分辨率
l = size(ts);
ts = imresize(ts, 1e4/l(1));
s = imbinarize(ts);%    二值图像
sumx = size(s);
sumy = sumx(2);
sumx = sumx(1);
[sdf1] = bwdist(s);   %   内部sdf
sdf1 = -sdf1;
[sdf2] = bwdist(~s);  %外部sdf
sdf = sdf2;
sdf(sdf2 == 0) = sdf1(sdf2 == 0);
end
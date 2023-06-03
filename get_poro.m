function [phi,phi_b] = get_poro(gca, gcf)
%   获取不包含死孔的孔隙度phi,包含死孔孔隙度phi_b
%   提供图窗gca和gcf,要求只有固体-孔隙两个颜色
set(gca,'position',[0,0,1,1]);
set(gcf,'position',[0,0,700,700]);
axis equal;
saveas(gca,'data\pic\get_poro_pic.bmp');
t = imread('data\pic\get_poro_pic.bmp');
ts = rgb2gray(t);
s = imbinarize(ts);%    二值图像
% figure;
% imshowpair(ts,s,'montage');% 画出二值图像，检查 
A = numel(s);
Asolid = A - sum(s,'all');   %   固体面积
phi_b = 1 - Asolid/numel(s);
st = imfill(~s,[1 1],4);    %   反转二值图像泛洪填充白色，于是只有死孔是黑色
% figure;
% imshow(st);
Adeadpore = A - sum(st,'all');
phi = 1 - (Adeadpore + Asolid)/A;

end
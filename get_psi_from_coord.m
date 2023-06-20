function [i, j, t] = get_psi_from_coord(x, y, Xlim, Ylim, sumx, sumy)
%   s是二值图像，固体为1，bulk为0，整张图像的范围为Xlim和Ylim,给定位置x,y,求出对应的像素点位置i，j
%   sumx,sumy为图像xy像素个数
%   t表示缩放倍率:实际一个像素的大小
if(x < Xlim(1) || x > Xlim(2))
    error('x给错了');
end
if(y < Ylim(1) || y > Ylim(2))
    error('y给错了');
end
sdx = (Xlim(2) - Xlim(1))/sumx;
j = floor((x-Xlim(1))/sdx) + 1;
sdy = (Ylim(2) - Ylim(1))/sumy;
i = floor((y-Ylim(1))/sdy) + 1;
if(i > sumy)
    i = i - 1;
end
if(j > sumx)
    j = j - 1;
end
i = sumy + 1 - i;
tx = (Xlim(2) - Xlim(1))/sumx;
ty = (Ylim(2) - Ylim(1))/sumy;
if(abs(tx - ty) > 1e-4)
    error('图片xy不对');
end
t = tx;
end
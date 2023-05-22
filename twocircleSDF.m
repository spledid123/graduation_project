function [sd, kx, ky, ispore] = twocircleSDF(x, x1, y1, x2, y2, r, h, len, flag)
%   对于六边形孔喉结构，将其分解为一段段直线，即两个圆来看待
%   对于两个连接在一起的圆，计算其SDF,点的位置为x = [x,y]
%   ispore==0表示在喉;ispore==1表示在孔
%   两个圆的圆心为[x1,y1],[x2,y2],半径为r,喉道的宽度为h,喉道长度为len
%   flag = 0表示横线，flag==1表示斜向右上,flag==-1表示斜向右下
%   首先归一到(0,0),(-r*(2+a),0)的两个圆（平移旋转）；再归一到(0,0)这个圆（对称）

% 归一化
t = r + len/2;%   六边形边长的一半
dx = x2;
dy = y2;
rx = [x(1) - dx; x(2) - dy];
l = 3^0.5*0.5;
if(flag == 1)
    rx = [0.5 l;-l 0.5]*rx;
elseif(flag == -1)
    rx = [0.5 -l;l 0.5]*rx;
end
kx = floor(x2/t);
ky = floor(y2/l/t/2);
if(rx(1)<-len/2-r)
    rx(1) = -(len + r * 2) - rx(1);
    kx = floor(x1/t);
    ky = floor(y1/l/t/2);
end

k = (r^2-(r-h/2)^2)/2/r/(r-h/2);
sr = -h/2 + k*r;%  小圆半径
sx = -r;
sy = k * sx;
th1 = pi-atan(k);
th2 = atan2(h/2,-r);
if(rx(1)==0&&rx(2)==0)
    phi = 0;
else
    phi = atan2(rx(2),rx(1));
end
ispore = 1;
if(phi<th1&&phi>-th1)
    sd = r - (rx(1)^2+rx(2)^2)^0.5;
elseif(phi<th2&&phi>0)
    sd = ((rx(1)-sx)^2+(rx(2)+sy)^2)^0.5-sr;
elseif(phi>-th2&&phi<0)
    sd = ((rx(1)-sx)^2+(rx(2)-sy)^2)^0.5-sr;
elseif(rx(1)>-r)
    sd = min(((rx(1)-sx)^2+(rx(2)+sy)^2)^0.5-sr, ((rx(1)-sx)^2+(rx(2)-sy)^2)^0.5-sr);
else
    ispore = 0;
    sd = min(h/2 - rx(2), rx(2) + h/2);
end
end
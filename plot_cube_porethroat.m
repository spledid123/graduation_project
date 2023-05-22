%   画出正方体排列的孔喉阵列
%   表示六边形的优势：粒子不会飞出去
hold on;
axis equal;
r = 1;
a = 0.25;%  喉道长度
R = 1.1;%   孔喉比
t = 2*r+2*a*r;
kx = -1:0;
ky = -1:0;
for j = 1:length(ky)
    for i=1:length(kx)
       plot_onecir_in_cube(t * kx(i), t * ky(j), r, 2 * (R - 1) * r, 2 * a * r);
    end
end


function [] = plot_onecir_in_cube(x, y, r, h, len)
%   画出正方形排列的孔喉阵列的一个圆
%   （x,y）表示圆心
%   r表示圆的半径
%   h为喉道宽度
%   len为喉道长度
    th1 = asin(h/2/r);
    th = th1:0.01:pi/2-th1;
    plot_pore_cir(x,y,r,th);
    x1 = r * cos(th1);
    x2 = r+len/2;
    plot([x+x1 x+x2],[y+h/2 y+h/2],'b');
    plot([x+h/2 x+h/2],[y+x1 y+x2],'b');
    th=th+pi/2;
    plot_pore_cir(x,y,r,th);
    plot([x-x1 x-x2],[y+h/2 y+h/2],'b');
    plot([x-h/2 x-h/2],[y+x1 y+x2],'b');
    th=th+pi/2;
    plot_pore_cir(x,y,r,th);
    plot([x-x1 x-x2],[y-h/2 y-h/2],'b');
    plot([x-h/2 x-h/2],[y-x1 y-x2],'b');
    th=th+pi/2;
    plot_pore_cir(x,y,r,th);
    plot([x+x1 x+x2],[y-h/2 y-h/2],'b');
    plot([x+h/2 x+h/2],[y-x1 y-x2],'b');
end
function []=plot_pore_cir(x,y,r,th)
%   画圆和临近孔，只画一部分，th给出圆的范围
    rx = x + r * cos(th);
    ry = y + r * sin(th);
    plot(rx,ry,'r');
end
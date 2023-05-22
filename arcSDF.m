function sd = arcSDF(x, c, th1, th2, r1, r2)
%   圆弧的SDF
%   圆弧的形状可以理解为圆环被两根角度为th1和th2的半径切断，且th1逆时针到th2之间的部分消失；
%   以这两条直边为直径加上半圆,保证半圆之间不会相碰
%   th的设置：x轴正方向为0，逆时针为正
%   x = [x, y]为点的位置,c = [cx, cy]为圆心,r1 < r2 为圆环半径
if(r1 > r2)
    error('圆弧半径设置错误');
end
th1 = th1 - floor(th1 / (2*pi)) * 2*pi;
th2 = th2 - floor(th2 / (2*pi)) * 2*pi;
if(th1 > th2)
    th1 = th1 - 2*pi;
end
%   旋转归一，令被切断的部分中点在x正半轴上,圆心在零点，圆弧关于x轴对称,点位于第一二象限
thmean = (th1 + th2)/2;
th = (th2 - th1)/2;
if(th < asin((r2 - r1)/(r1 + r2)))
    error('圆弧角度设置错误，半圆相碰');
end
r = x - c;
r = [cos(thmean), sin(thmean); -sin(thmean), cos(thmean)] * r';
r = r';
r(2) = abs(r(2));
%   圆弧部分由th-pi;o-th是空白的
thr = atan2(r(2), r(1));
if(thr > th)
    t1 = circleSDF(r, [0,0], r1);
    t2 = circleSDF(r, [0,0], r2);
    if(t2 > 0)
        sd = t2;
    elseif(t1 < 0)
        sd = -t1;
    else
        sd = -min(abs(t1),abs(t2));
    end
else
    sd = circleSDF(r, (r1 + r2)/2*[cos(th),sin(th)], (r2 - r1)/2);
end
end
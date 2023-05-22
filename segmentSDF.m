function sd = segmentSDF(x, a, b)
%   计算线段AB的SDF
%   x = [x, y]为点X的坐标,a = [ax, ay]；b = [bx, by]为线段的两个端点A,B
v = x - a;% 相对于端点A的坐标,记为AX
u = b - a;% 线段矢量AB
t = max(min(dot(v,u)/dot(u,u), 1), 0);
d = v - u * t;
sd = norm(d);
%   将空间分为三个部分
%   假设A为零点，AB的方向为正x轴
%   若点在第二三象限x < 0, 则t = 0, d = v = AX
%   若x > 0,且x < |AB|,且AX与AB的夹角为theta ,则t = cos(theta),于是d为垂线
%   若x > |AB|, 则t = 1, d = BX

end
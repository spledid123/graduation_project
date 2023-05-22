function sd = triangleSDF(x, a, b, c)
%   求三角形ABC的SDF
%   x = [x, y]为点X的坐标；a = [ax, ay], b = [bx, by], c = [cx, cy]为A，B，C的坐标
d = min([segmentSDF(x, a, b), segmentSDF(x, a, c), segmentSDF(x, b, c)]);%  X到三角形的最短距离
t = sign(dim2cross(b - a, x - a)) == sign(dim2cross(c - b, x - b)) &&...
    sign(dim2cross(b - a, x - a)) == sign(dim2cross(a - c, x - c));
%   二维向量的叉乘w×v为在z轴上的向量
%   叉乘的正负能够反映两条向量的sin(theta),于是若AB×AC与AB×AD同方向(z的正向/负向)，则CD在线段AB的同侧
%   若在三角形ABC内部，则必然在AB,BC,CA的同侧；反之则不然
%   于是，若在ABC内部，则三个叉乘要不都大于0，要不都小于0
if(t)
    sd = -d;
else
    sd = d;
end

end

function s = dim2cross(a, b)
%   二维向量的叉乘a×b为在z轴上的向量
%   cross为三维向量的叉乘
%   s为a×b的z轴分量,a,b为二维向量
ra = a;
ra(3) = 0;
rb = b;
rb(3) = 0;
s = cross(ra, rb);
s = s(3);
end
function sd = circleSDF(x,c,r)
%   计算到某个圆的距离
%   x=[x,y]为点的位置
%   c=[cx,cy]为圆心位置
%   r为半径 

ux=x(1)-c(1);
uy=x(2)-c(2);
sd = sqrt(ux*ux+uy*uy)-r;
end
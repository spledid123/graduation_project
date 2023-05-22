function n = gradient(x , R)
%   求边界的法向量n=(nx,ny)
%   x=[x,y]为边界点的坐标，R为形状参数，为函数scene的参数
EPSILON = 1e-6;

n(1) = (scene([x(1) + EPSILON x(2)],R) - scene([x(1) - EPSILON x(2)],R)) * 0.5 / EPSILON;
n(2) = (scene([x(1) x(2) + EPSILON],R) - scene([x(1) x(2) - EPSILON],R)) * 0.5 / EPSILON;

len = sqrt(n(1)^2+n(2)^2);
n = n ./ len;
end
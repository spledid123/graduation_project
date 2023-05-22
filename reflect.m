function r = reflect(i , n)
%   反射，改变方向
%   i=[ix,iy]为入射法向量
%   n=[nx,ny]为表面法向量
%   r=[rx,ry]为出射法向量

fcosin = @(x) asin(2 * x -1);
%   余弦分布,与入射角无关
alpha = fcosin(rand());
r = n * [cos(alpha) sin(alpha) ; -sin(alpha) cos(alpha)];

end
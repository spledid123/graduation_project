function n = gradient(x, R, EPSILON_SD, varargin)
%   求边界的法向量n=(nx,ny)
%   x=[x,y]为边界点的坐标，R为形状参数，为函数scene的参数
p = inputParser;            % 函数的输入解析器
addOptional(p,'pm',0);
parse(p,varargin{:});
ppm = p.Results.pm;


n(1) = (scene([x(1) + EPSILON_SD x(2)],R,ppm) - scene([x(1) - EPSILON_SD x(2)],R,ppm)) * 0.5 / EPSILON_SD;
n(2) = (scene([x(1) x(2) + EPSILON_SD],R,ppm) - scene([x(1) x(2) - EPSILON_SD],R,ppm)) * 0.5 / EPSILON_SD;

len = sqrt(n(1)^2+n(2)^2);
n = n ./ len;
end
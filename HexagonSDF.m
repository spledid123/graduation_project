function sd = HexagonSDF(x, c, r)
%   正六边形的sdf
%   x表示位置,c表示中心点,r表示中心点到顶点的距离,th表示逆时针旋转距离,n表示n边形
r = r/2*3^0.5;
p = x - c;
k = [-cos(pi/6) sin(pi/6)];
p = abs(p);
p = p - 2*min(dot(k, p), 0)*k;
p = p - [clamp(p(1), -tan(pi/6)*r, tan(pi/6)*r) r];
sd = norm(p) * sign(p(2));
end
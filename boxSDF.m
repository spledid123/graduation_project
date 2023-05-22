function sd = boxSDF(x, c, s, th)
%   计算点到矩形的SDF
%   x = [x, y]为点的坐标;c = [cx,cy]为矩形中心坐标，s = [sx,sy]为矩形的半长轴与半宽轴
%   th = 0时长轴为x轴，宽轴为y轴；th表示矩阵逆时针转的角度
%   考虑到计算误差，为矩形做倒角，默认四个顶点有90°小圆弧，半径为Delta
Delta = 1e-3;
r = [cos(th), sin(th); -sin(th), cos(th)] * (x - c)';%   相对于th = 0矩形中心点的坐标
r = r';
r = abs(r);%    都到第一象限
d = r - s; %   相对于顶点的坐标,d = [dx, dy]
ax = max(0, d(1));%   如果没有超过某条轴，记为0
ay = max(0, d(2));
sd = min(max(d(1), d(2)), 0) + (ax ^ 2 + ay ^ 2) ^ 0.5;
%   看第一象限，
%   若在矩形外部
%   如果dx > 0, dy < 0,max(dx, dy) = dx > 0,min(max(dx, dy), 0) = 0; ax = dx,ay = 0 -> sd = dx；
%   若dx < 0, dy > 0,同理,sd = dy
%   若dx > 0, dy > 0,则min(max(dx, dy), 0) = 0；sd = (dx ^ 2 + dy ^ 2) ^ 0.5;
%   若在边上，dx = 0,dy < 0,则min(max(dx, dy), 0) = 0，ax = 0, ay = 0, sd = 0；在另一条边上同理
%   若在内部,dx < 0, dy < 0,则min(max(dx, dy), 0) = max(dx,dy), ax = 0, ay = 0, sd = max(dx,dy)<0;
if((d(1) > -Delta) && (d(2) > -Delta))
    sd = circleSDF(r, s - Delta, Delta);
end
%   进行倒角，(-Delta < dx) && (-Delta < dy) 的SDF计算都将用圆的方案计算
%   若用第一象限的坐标r计算,圆心为(sx - Delta, sy - Delta), 半径为Delta
end

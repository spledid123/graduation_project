function phi = get_poro(a, color)
%   获取孔隙度phi
%   提供图窗gca,要求只有固体-孔隙两个颜色
%   color为填充固体的颜色，三元组
frame = getframe(a);
f = frame.cdata;
su = size(f);
sumi = su(1);
sumj = su(2);
phi = 0;
for i = 1:sumi
    for j = 1:sumj
        if(f(i,j,1) == color(1) && f(i,j,2) == color(2) && f(i,j,3) == color(3))
            phi = phi + 1;
        end
    end
end
phi = 1 - phi/(sumi*sumj);
end
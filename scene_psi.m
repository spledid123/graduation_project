function sd = scene_psi(x, y, sdf, Xlim, Ylim, sumx, sumy)
%   通过像素sdf得到实际sdf
[i, j, t] = get_psi_from_coord(x, y, Xlim, Ylim, sumx, sumy);
if(i < 1 || i > sumx || j < 1 || j > sumy)
    error('xy范围不对,x=%f,y=%f,i=%d,j=%d',x,y,i,j);
end
sd = sdf(i,j)*t;
end
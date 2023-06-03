function sd = scene_psi(x, y, sdf, Xlim, Ylim, sumx, sumy)
%   通过像素sdf得到实际sdf
[i, j, t] = get_psi_from_coord(x, y, Xlim, Ylim, sumx, sumy);
sd = sdf(i,j)*t;
end
function sd = get_poro_six(R, a)
%   给定六边形孔喉结构的参数，返回孔隙度
sd = (2*pi+12*a*(R-1))/(8*3^0.5*(1+a)^2);
end
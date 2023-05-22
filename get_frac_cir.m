%   满足分形规律的圆阵列
kx = -10:9;
ky = -10:9;
numx=length(kx);
numy=length(ky);
R = [1.1,1.1];
axis off
axis equal 
[NR,Df] = get_fractal_cir_R(0.1,2.2,400,1.5,"NUM");
for i=1:numx
    x0(i)=R(1)*(2*kx(i)+1);
end
for i=1:numy
    y0(i)=R(2)*(2*ky(i)+1);
end
hold on;
for i=1:numx
    for j=1:numy
        plot_cir(x0(i),y0(j),NR((i-1)*numx+j));
    end
end

function [R, D] = get_fractal_cir_R(Rmin, Rmax, NUM, Df, flag)
%   生成满足分形分布的圆半径
%   Rmin,Rmax分别为半径的最小/最大值
%   P[X>=R] = (Rmax/R)^Df/NUM
%   数量有限,总数NUM = (Rmax/Rmin)^Df;
%   Df为分形维数，2维小于2；3维小于3
%   R为数量为NUM的数组
%   若flag=='NUM',则输入NUM,输出Df;若flag=='Df',则输入Df,输出NUM

if(flag=="Df")
NUM = (Rmax/Rmin)^Df;
D = NUM;
elseif(flag == "NUM")
Df = log(NUM)/log(Rmax/Rmin);
D = Df;
end
i = 1;
while(i < NUM + 1)
    a = rand();
    R(i) = Rmax/((a*NUM)^(1/Df));
    if(R(i)>Rmax||R(i)<Rmin)
        i = i - 1;
    end
    i = i + 1;
end


end

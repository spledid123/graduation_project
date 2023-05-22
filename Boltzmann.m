function s = Boltzmann(T)
%   输出出射的速度大小，T为温度
%%  3D玻尔兹曼分布,2D泄流分布
%   产生满足玻尔兹曼分布的随机数
%   给定绝对温度,pi为速度
% 
k = 1.380649e-23;
NA = 6.02e23;
m = 18e-3/NA;

a = sqrt(k * T / m);
norm_cdf = @(x) erf(x / sqrt(2) / a) - x * exp(-x * x / (2 * a * a)) / a * sqrt(2 / pi);
%   分布的cdf

p = rand();
s = 150;
if (norm_cdf(s) < p)
    while (norm_cdf(s) < p)
        s = s + 0.01;
    end
end
if (norm_cdf(s) > p)
    while (norm_cdf(s) > p)
        s = s - 0.1;
    end
end
%%  2D玻尔兹曼分布
%   f(v)=2Avexp(-Av^2);F(v)=1-exp(-Av^2);v=sqrt(ln(1-F(v))/(-A)),F(v)=rand();
% k = 1.380649e-23;
% NA = 6.02e23;
% m = 18e-3/NA;
% 
% A = m/2/k/T;
% p = rand();
% s = ((log(1-p))/-A)^0.5;

%%  匀速
% s=100;


end
function fz = fun_exp_bol_pdf(z)
%   概率密度函数
%   指数分布/玻尔兹曼分布的和
%   f(z)=f(X+tY),X~Bolt;Y~exp
lam1 = 0.1;
lam2 = 0.03;
A = 0.01;
% v = 1000;
t = 10;
k = 1.380649e-23;
NA = 6.02e23;
m = 18e-3/NA;
% T = v^2*pi*m/8/k;
T=40;
B = m / 2 / k / T;
A = 4*B^1.5/pi^0.5;

% 指数+玻尔兹曼
% fun=@(x)exp(-(z-x)./lam-B./(x./t).^2)./(x./t).^4.*A./lam /t;
% % fz = quadgk(fun,0,z);
%   指数分布
% fun=@(x)exp(-x/lam)/lam;
% fz=fun(z);
%   指数pdf+玻尔兹曼pdf
fun = @(x)exp(-x/lam1)/lam1*A+exp(-x/lam2)/lam2*(1-A);
fz=fun(z);
end
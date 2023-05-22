function s = T_absorb()
%   返回与边界碰撞后吸附的时间
%%  不考虑吸附
s = 0.00;
%%  指数分布的吸附时间
% T = 0;% 平均吸附时间
% x = rand();
% fun = @(x)-log(x)/T;
% s = fun(x);

%%  玻尔兹曼分布的吸附时间（不考虑该分布）
% k = 1.380649e-23;
% NA = 6.02e23;
% m = 18e-3/NA;
% T = 0.0008504575;%  平均时间为1
% 
% a = sqrt(k * T / m);
% norm_cdf = @(x) erf(x / sqrt(2) / a) - x * exp(-x * x / (2 * a * a)) / a * sqrt(2 / pi);
% %   分布的cdf
% 
% p = rand();
% s = 1;
% if (norm_cdf(s) < p)
%     while (norm_cdf(s) < p)
%         s = s + 1e-4;
%     end
% end
% if (norm_cdf(s) > p)
% 
%     while (norm_cdf(s) > p)
%         s = s - 1e-4;
%     end
% end


end
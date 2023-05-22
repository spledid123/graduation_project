%   用有限差分的办法直接求解一个简单的1D扩散方程
%   Δc/Δt=D(Δ^2c)/(Δx)^2,x∈[0 L],0处是封闭边界,L处是封闭或定浓度边界(模拟吸附)
%   初始条件:c(x)=2c0 - (2c0/L)*x
%   差分：c_new(n)-c(n)=D*dt/(dx)^2*(c(n+1)+c(n-1)-2*c(n))
%   结论：吸附会显著降低扩散的平衡时间，封闭边界平衡时间大概为L^2/2/D;定浓边界大概是2L^2/D
L = 1e-4;
D = 73.6e-5;
c0 = 1e-2;

n = 100;    % 差分步数
step = 20e4; %   时间长度或迭代次数
c = zeros(1, n+1);%    c(x)
c_0 = 0;    %   左端边界
c_L = 0;    %   右端边界
dt = L^2/D/1e5;
dx = L/n;
c_new = zeros(1, n+1);
c_sta = [];%    初始条件
c_r = 0;    %   c(n+1)
c_l = 0;    %   c(n-1)
cc = [];    %   临时记录某时刻的浓度，方便画图

%   初始化
f = @(x)2*c0-(2*c0/L)*x;
for i = 1:n+1
    c(i) = f((i-1)*dx);
end
c_sta = c;
c_0 = c(1);
c_L = c(n+1);
%   求解
for i = 1:step % 迭代次数
    for j = 1:n+1   %   求解
        if(j == 1)
            c_l = c_0;
            c_r = c(2);
        elseif(j == n+1)
            c_l = c(n);
            c_r = c_L;
        else
            c_l = c(j - 1);
            c_r = c(j + 1);
        end
        c_new(j) = c(j) + D * dt / (dx)^2 * (c_r + c_l - 2 * c(j));
    end
    %   边界条件
    %   0处封闭边界
    c_0 = c(1);
    %   L处
    c_L = c(n + 1);     %   封闭条件
    c_L = 0;            %   定浓度
    %   赋值
    c = c_new;
    if(i == 5e4)
        cc = c;
    end
end
figure;
hold on;
box on;
plot((0:n)*dx/L, c_sta/c0,'.b', 'displayname','初始时刻');
t = 1e4;
t = num2str(t);
t = [t '时刻' '右端封闭'];
plot((0:n)*dx/L, cc/c0,'-k', 'displayname',t);
t = step;
t = num2str(t);
t = [t '时刻' '右端封闭'];
plot((0:n)*dx/L, c/c0,'--k', 'displayname',t);
title('浓度-位置','fontsize',14);
xlabel('位置','fontsize',14);
ylabel('浓度','FontSize',14);




%   多粒子1维随机游走，考虑吸附
%   画出粒子在每一步的吸附概率与扩散系数的关系
%   线性关系
clc;
clear;
pn=0:0.1:1;% 吸附概率
n = length(pn);
t=0;
for k=1:n % 每个概率进行模拟
    p=pn(k);%   概率
    N_T = 1000;%步数
    N = 500;%粒子数
    s2N = zeros(1,N_T+1);%  每一步的均方位移

    for i = 1:N %   多粒子
        s = 0; %    粒子的位置
        for j = 1:N_T % 游走
            m = rand();
            if(m<(1-p)/2)
                t=1;
            elseif(m<1-p)
                t=-1;
            else
                t=0;
            end
            s=s+t;
            s2N(j+1)=s2N(j+1)+s^2;
        end
    end
    s2N = s2N/N;
    D(k)=(s2N(N_T+1)-s2N(1))/N_T;%  扩散系数
end
plot(pn,D);%    画出吸附概率与扩散系数，线性    


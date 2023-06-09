%   以及get_porethroat_para.m得到的参数t0,pp,pt,Tp,Tt
%   根据孔喉双重结构随机游走模型计算六边形孔隙网络的扩散系数
R = 1.1;
a = 0.25;
Tp = 1.4334e-7;             %   在孔内的单步时间
Tt = 2.0702e-8;              %   在喉内的单步时间
% Tp = 1e-5;% 平均吸附时间
% Tt = 1e-5;
pp = 0.0957;            %   单次碰撞孔离开孔的概率
pt = 0.323;               %   单次碰撞喉离开喉的概率
fp = -log(1-pp)/Tp;
ft = -log(1-pt)/Tt;
T0 = 1.056e-7;  %   单步时间
% T0 = 1e-5;
Pp = 1-exp(-fp*T0); %   单步概率
Pt = 1-exp(-ft*T0);
%   x方向扩散系数等于方差除以单位时间
a = (Pp*Pt) / (2*Pt + 3*Pp);
si = (2*5*a/3 + 8*5*a/6);
dr = (1 + a)*1e-5/2;
si2 = si * dr^2;%   方差
Dx = si2 / 4 / T0;% 扩散系数,2D除以4
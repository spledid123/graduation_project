%   画速度云图，可能会有归一化
%   画需要的多孔介质

filenamesta = 'data\diffusion_R_1e-5_1.1_dT_1e-5\rxT_circle_T_';
filenamemid = num2str(0.00041);
filenameend = '_.txt';
filename = strcat(filenamesta,filenamemid,filenameend);
A = readtable(filename);
x = A.rx; % 提取x坐标
y = A.ry; % 提取y坐标
v = A.v; % 提取速度大小

% 可能修正半径
% x = x * 1e5;y = y * 1e5;
%   圆阵列位置归一化
Rp = 1.1;
kx = floor((x / Rp - 1) / 2);
ky = floor((y / Rp - 1) / 2);
X = x - 2 * Rp * (kx + 1);
Y = y - 2 * Rp * (ky + 1);
%   开放空间圆阵列归一化
% r = 500;
% Rp = 1.1;
% nump = 30;
% t = 2 * Rp;
% tx = t;
% ty = t;
% kx = floor(x/tx);
% ky = floor(y/ty);
% x = x - kx * tx;
% y = y - ky * ty;
%   筛选粒子
% t_pore = (x < Rp * nump) & (x > -Rp * nump) & (y < Rp * nump) & (y > -Rp * nump);
% % t_pore = [];
% % for k = 1:length(v)
% %     t_pore(k) = triangleSDF([x(k) y(k)], [-115*3^0.5,-115], [115*3^0.5,-115], [0, 230]) < 0;
% % end
% t_cir = (x.^2 + y.^2).^0.5 < 60-1e-3;
% j = 1;
% for i = 1:length(v)
%     if(t_cir(i) == 1)
%         V(j) = v(i);
%         X(j) = x(i);
%         Y(j) = y(i);
%         j=j+1;
%     end
% end
V = v;

% 计算不同位置的粒子平均速度
[x,y] = meshgrid(min(X):0.02:max(X),min(Y):0.02:max(Y)); % 生成方格内的网格点
v = griddata(X,Y,V,x,y); % 使用griddata函数插值得到每个网格点的速度大小
% 绘制云图，使用颜色表示速度大小
pcolor(x,y,v); % 使用pcolor函数绘制云图
shading interp; % 平滑颜色过渡
colorbar; % 显示颜色条
axis equal;
hold on;
xlabel('x'); % 标注x轴
ylabel('y'); % 标注y轴
title('Average velocity of particles at different positions'); % 标注标题
xlim([-Tlim Tlim]);
ylim([-Tlim Tlim]);

%   画多孔介质
Tlim = 1.1;
plot_cir(Tlim,Tlim,1);
plot_cir(Tlim,-Tlim,1);
plot_cir(-Tlim,Tlim,1);
plot_cir(-Tlim,-Tlim,1);

% plot_cir(0,0,500);
% numr = 30;    % 右边圆阵列圆数量
% numl = 50;    % 左边圆阵列圆数量
% Rr = 1;       % 右边圆阵列圆半径
% Rl = 1.5;      % 左边圆阵列圆半径
% tr = 2;      % 右边圆阵列圆半径+喉道宽度一半
% tl = 2.2;      % 左边圆阵列圆半径+喉道宽度一半
% tsdr = -60;      %  圆阵列所在正方形与区块的间隔
% tsdl = 0;
% %   右边区块
% kx = 0:numr - 1;
% ky = 0:numr - 1;
% for i=1:numr
%     x0(i)=tr*(2*kx(i)+1) + tsdr;
% end
% for i=1:numr
%     y0(i)=tr*(2*ky(i)+1) + tsdr;
% end
% for i=1:numr
%     for j=1:numr
%         plot_cir(x0(i),y0(j),Rr);
%     end
% end
% xlim([-500 500]);
% ylim([-500 500]);

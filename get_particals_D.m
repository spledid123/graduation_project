%  根据多粒子轨迹，画<R^2>-t，计算扩散系数
%   D = <R^2>/t
%   线性回归
filenamesta = "data\diffusion_R_1e-5_1.1_dT_1e-5\rxT_circle_T_";
filenameend = "_.txt";
N_T = 1e-5;
ii = 1:43;
r2x = zeros(1,1e3);
r2y = zeros(1,1e3);
r2R = zeros(1,1e3);
th = pi / 4; %    计算某个特定方向角的均方位移
r2thx = zeros(1,1e3);
r2thy = zeros(1,1e3);
tic;
for i = ii
    filename = strcat(filenamesta,num2str(i * N_T),filenameend);
    A = readtable(filename);
    x = A.rx;
    y = A.ry;
    len = length(x);
    r2x(i) = sum(x.^2) / len;% x方向均方位移
    r2y(i) = sum(y.^2) / len;% y方向均方位移
    r2R(i) = r2x(i) + r2y(i);      % 均方位移
    A = [x y] * [cos(th) -sin(th);sin(th) cos(th)];
    x = A(:,1);
    y = A(:,2);
    r2thx(i) = sum(x.^2) / len;%   th方向角的均方位移
    r2thy(i) = sum(y.^2) / len;%   th方向角垂直方向的均方位移
end
toc;
figure;
hold on;
box on;
xlabel('时间/s','fontsize',14);
ylabel('均方位移/m^2','fontsize',14);
title('均方位移-时间','fontsize',14);
plot(ii*N_T,r2R(ii),'DisplayName','均方位移');
plot(ii*N_T,r2x(ii),'DisplayName','x方向');
plot(ii*N_T,r2y(ii),'DisplayName','y方向');
plot(ii*N_T,r2thx(ii),'DisplayName','th方向');
plot(ii*N_T,r2thy(ii),'DisplayName','th垂直方向');
px = polyfit(ii*N_T,r2x(ii),1);
py = polyfit(ii*N_T,r2y(ii),1);
pthx = polyfit(ii*N_T,r2thx(ii),1);
pthy = polyfit(ii*N_T,r2thy(ii),1);
p = polyfit(ii*N_T,r2R(ii),1);
fprintf("px = %f, py = %f, p = %f\npthx = %f, pthy = %f",px(1), py(1), p(1), pthx(1), pthy(1));
%%   画出扩散系数和R的关系
D1 = [0.006446 0.064351 0.227403 1.550708];
D1 = D1/1e-4;
D2 = [0.000736 0.007454 0.018093 0.144077];
D2 = D2/1e-5;
R = [1.1 1.5 2 5];
x = (2^0.5.*R - 1)./(R - 1);%   圆阵列的实际孔喉比
figure;
hold on;
box on;
set(gca,'xscale','log')
yyaxis left;
plot(D1,R,'o-b','displayname','r=1e-4m');
plot(D2,R,'o-r','displayname','r=1e-5m');
xlabel('扩散系数/半径','fontsize',14);
ylabel('R——(圆直径+喉道宽度)/圆直径','fontsize',14);
title('修正扩散系数-R','fontsize',14);
yyaxis right;
plot(D1,x,'o--b','displayname','r=1e-4m');
plot(D2,x,'o--r','displayname','r=1e-5m');
ylabel('孔喉比','FontSize',14);
ax = gca;
yyaxis right;
ax.YDir = 'reverse';
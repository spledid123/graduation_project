%   运用并行parfor得到一系列粒子的路径库，需要调整scene以改变形状，调整reflect调整反射模式
%   需要给出每个粒子的形状参数
x0 = [0 0];%    初始位置
gx0 = [0 0];%   初始速度方向
R = [];%     形状系数

s = [];%    经过一定次数的碰撞后，粒子的位移大小
x = [];%    经过一定次数的碰撞后，粒子的x位移
y = [];%    经过一定次数的碰撞后，粒子的y位移
gx = [];

num_par = 500;% 粒子数 
N_pmax = 2e3;%  最大碰撞次数
%   不同粒子的形状参数
for i = 1:num_par
    R(i,:) = [1.1 1.1];
end

%   文件名
filenamesta = 'data\path_2000_six(1)\rxT_circle_onepartical_sixcir_';
filenameend = '_.txt';
file1 = 'R_';
file2 = '_a_';

parfor m = 1:num_par
    filenamemid = int2str(m);
    filename = strcat(filenamesta,file1,num2str(R(m,1)),file2,num2str(R(m,2)),'_',filenamemid,filenameend);
    th = rand()*2*pi;
    gx0 = [cos(th) sin(th)];
    [filename, rx, gx, s(m)] = fun_partical_one(x0, gx0, N_pmax, filename, R(m,:));
    x(m) = rx(1);
    y(m) = rx(2);
end
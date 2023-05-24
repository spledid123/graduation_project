%   画出粒子不同时刻的位置，并生成gif
%   该程序的逻辑是：gif的生成来自frame生成的img,而我们可以从画布gcf(也就是程序运行时的显示)生成frame.
%   我这里的frame其实是两个画布合并而来，做了一些处理——也就是说，显示的画布和gif其实不同。

clear;
FileName = 'data\pic\space.gif';% 创建文件名
axis off;
axis equal;
hold on;
set(gcf,'color','w');%  画布背景
%   图的范围
R = 500;
xlim([-R,R]);
ylim([-R,R]);

%   画多孔介质的图，可以从plot_partical_path.m搬运
% plot([200 200 -200 -200 200],[200 -200 -200 200 200],'r--','LineWidth',2);
%  x = 0左右两边的区块,圆阵列
% numr = 100;    % 右边圆阵列圆数量
% numl = 0;    % 左边圆阵列圆数量
% Rr = 1;      % 右边圆阵列圆半径
% Rl = 1;      % 左边圆阵列圆半径
% tr = 2;      % 右边圆阵列圆半径+喉道宽度一半
% tl = 2;      % 左边圆阵列圆半径+喉道宽度一半
% tsdr = -200;      %   区块间的间隔
% tsdl = 0;
% kx = 0:numr - 1;%   右边
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
% kx = 0:numl - 1;%   左边
% ky = 0:numl - 1;
% for i=1:numl
%     x0(i) = - (tl*(2*kx(i)+1)) - tsdl;
% end
% for i=1:numl
%     y0(i) = tl*(2*ky(i)+1) + tsdl;
% end
% for i=1:numl
%     for j=1:numl
%         plot_cir(x0(i),y0(j),Rl);
%     end
% end
%   圆随机场
% c = load('data\cir.mat');
% c = c.cir;
% for i = 1:1300
%     plot_cir(c(i,1),c(i,2),c(i,3));
% end
% plot([-106 * 3^0.5, 106 * 3^0.5],[-106 -106],'r');
% plot([-106 * 3^0.5, 0],[-106 212],'r');
% plot([0, 106 * 3^0.5],[212 -105],'r');


%   多孔介质的frame
frame_cir = getframe(gcf);
cla();% 清除画布(gcf)图像
N_T = 1;%   delta T
% 循环遍历每个表格
for i = 1:1:30
    % 从表格中提取粒子的x和y坐标
    filenamesta = 'data\bulk_R_500_N_200000_dT_2000\rxT_circle_T_';
    filenamemid = num2str(N_T * i);
    filenameend = '_.txt';
    filename = strcat(filenamesta,filenamemid,filenameend);
    A = readtable(filename);
    x = A.rx;
    y = A.ry;
    %   筛选粒子
    t = (x.^2 + y.^2 - R ^ 2) < -1e-3;
    x = x .* t;
    y = y .* t;
    % 绘制粒子图
    plot(x,y,'b.');
    plot(0,0,'w.');
    box on;
    title(['Time = ',num2str(i * N_T)])
    % 获取当前帧
    frame = getframe(gcf);
    % 将帧转换为图像
    frame.cdata = min(frame.cdata,frame_cir.cdata);%    保留圆的像素
    img = frame2im(frame);
    % 将图像转换为索引颜色
    [imgind,cm] = rgb2ind(img,256);
    cla();% 清除画布图像
    % 写入GIF文件
    if i == 1
        % 对于第一帧，开始“LoopCount”
        imwrite(imgind,cm,FileName,'gif','LoopCount',Inf,'DelayTime',0.01);
    else
        % 对于其他帧，追加到文件中
        imwrite(imgind,cm,FileName,'gif','WriteMode','append','DelayTime',0.01);
    end
end
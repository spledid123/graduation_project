Y=[];YY=[];NY=[];NYY=[];GY=[];GYY=[];GNY=[];GNYY=[];

a=gca;a=a.Children;
y=a(1).YData(3:37);yy=a(2).YData(3:37);
a=gca;a=a.Children;
ny=a(1).YData(3:37);nyy=a(2).YData(3:37);

y=y';yy=yy';ny=ny';nyy=nyy';
str = '正方形阵列(R=1.5)';
str2 = [str '内'];
str1 = [str '外'];
g = repmat({str1}, length(y), 1); %   外粒子数
GY = [GY;g];
% GY = g;
g = repmat({str1}, length(ny), 1); %  外数密度
GNY = [GNY;g];
% GNY=g;
g = repmat({str2}, length(yy), 1); %  内粒子数
GYY = [GYY;g];
% GY=g;
g = repmat({str2}, length(nyy), 1); % 内数密度
GNYY = [GNYY;g];
% GNYY=g;
Y = [Y;y];
YY = [YY;yy];
NYY = [NYY;nyy];
NY = [NY;ny];

figure;
ax1 = gca;
boxplot(ax1, Y, GY);
title('外粒子数');
figure;
ax2 = gca;
boxplot(ax2, YY, GYY);
title('内粒子数');
figure;
ax3 = gca;
boxplot(ax3, NY, GNY);
title('外数密度');
figure;
ax4 = gca;
boxplot(ax4, NYY, GNYY);
title('内数密度');


% boxplot(ax, y);
% hold on;
% group = repmat(str1, size(ny,1), 1); 
% boxplot(ax1, ny);
% hold on;
% str2 = [str '内'];
% group = repmat(str2, size(yy,1), 1); 
% boxplot(ax, y);
% group = repmat(str2, size(nyy,1), 1); 
% boxplot(ax1, ny);

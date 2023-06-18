%   配套analyse_particals.m中关于开放空间内多孔介质内外粒子数+数密度的分析
%   将以下语句分别贴到命令行窗口，把数密度信息存入matlab.mat
%   画箱型图

Y=[];YY=[];NY=[];NYY=[];GY=[];GYY=[];GNY=[];GNYY=[];

%   达到稳定的范围
%%
rr = 3:42;
%   粒子数
a = gca; a = a.Children;
y = a(1).YData(rr);yy = a(2).YData(rr);
%%   数密度
a = gca;a = a.Children;
ny = a(1).YData(rr);nyy = a(2).YData(rr);
%%
load('matlab.mat');
y = y';yy = yy';ny = ny';nyy = nyy';
str = 'pm18';
str2 = [str '内'];
str1 = [str '外'];
NYname{length(NYname) + 1} = struct('name', str1, 'sta', length(GNY) + 1, 'end', length(GNY) + length(ny));
NYYname{length(NYYname) + 1} = struct('name', str2, 'sta', length(GNYY) + 1, 'end', length(GNYY) + length(nyy));
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

%%   画箱型图
figure;
ax1 = gca;
boxplot(ax1, Y, GY);
title('外粒子数');
figure;
ax2 = gca;
boxplot(ax2, YY, GYY);
title('内粒子数');
%%  画内外数密度箱型图
figure;
subplot(1,2,1);
ax3 = gca;
boxplot(ax3, NY, GNY);
title('外数密度');
ylim([0 0.018]);
subplot(1,2,2);
ax4 = gca;
boxplot(ax4, NYY, GNYY);
title('内数密度');
ylim([0 0.018]);
%%
save('matlab.mat','GNYY','NYY',"GNY","NY","GYY","GY","YY","Y","NYname","NYYname");
%%  画指定的数密度箱型图
load('matlab.mat');
i = [13 14];
nyyi = [];
for j = 1:length(i)
    nyyi = [nyyi NYYname{i(j)}.sta:NYYname{i(j)}.end];
end
figure;
ylim([0 0.018]);
ax5 = gca;
boxplot(ax5, NYY(nyyi), GNYY(nyyi));
title('内数密度');
ylim([0 0.018]);
function [xxfit,yy] = plot_distribution(kTTT,varargin)
%   给出数组kTTT的概率分布函数
%   默认连续,第二个可选参数为flag,flag=0为连续，flag=1为离散整数，连续的纵坐标为频率密度,离散整数时为频率
%   若为连续，则区间取edges个,属于名称-数值，默认100;
%   实例：plot_distribution(kTTT);plot_distribution(kTTT,0,'edges',360);

p = inputParser;            % 函数的输入解析器
fcn = @(x)length(x) > 1 && isnumeric(x);
ffcn = @(x) x == 1 || x ==0;
addRequired(p,'width',fcn);
addOptional(p,'flag',0,ffcn);%    默认连续量
addParameter(p,'edges',100);      % 默认100个区间
parse(p,kTTT,varargin{:});       % 对输入变量进行解析，如果检测到前面的变量被赋值，则更新变量取值
flag = p.Results.flag;
edges = p.Results.edges + 1;
if(flag == 0)
    figure;
    x=linspace(min(kTTT),max(kTTT),edges); %计算各区间中点
    yy=histogram(kTTT,x,'Normalization','pdf'); %计算各个区间的个数
    xx=yy.BinEdges;
    yy=yy.Values;
    for i=1:length(xx)-1
        xx(i)=(xx(i)+xx(i+1))/2;
    end
    hold on;
    xxfit = xx(1:length(xx)-1);
    plot(xxfit,yy);
else
    xxfit=1:max(kTTT);
    yy = zeros(1,max(kTTT));
    for i=1:length(kTTT)
        yy(kTTT(i)) = yy(kTTT(i)) + 1;
    end
    yy = yy/sum(yy);
end
end
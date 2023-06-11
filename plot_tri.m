function [] = plot_tri(a, b, c, varargin)
%   画三角形
%   a,b,c是三个顶点
%   flag为可选参数,flag == 0默认,flag == 1填充
%   名称-数值对:var为plot的参数；fillvar为填充颜色
p = inputParser;            % 函数的输入解析器
fcn = @(x) x == 0 || x == 1;
addOptional(p,'flag',0,fcn);
addParameter(p,'var',"r-");   
addParameter(p,'fillvar',"r");      % 默认参数
parse(p,varargin{:});
if(p.Results.flag)
    fill([a(1) b(1) c(1)], [a(2) b(2) c(2)], p.Results.fillvar,'LineStyle','none');
else
    plot([a(1) b(1) c(1) a(1)], [a(2) b(2) c(2) a(2)],p.Results.var);
end
end
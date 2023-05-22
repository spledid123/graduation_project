function [] = plot_rec(c,s,th,varargin)
%   画长方形
%   c=[cx,cy]为中心，s=[sx,sy]为半长半宽,th为逆时针旋转角度
%   flag为可选参数,flag == 0默认,flag == 1填充
%   名称-数值对:var为plot的参数；fillvar为填充颜色
p = inputParser;            % 函数的输入解析器
fcn = @(x) x == 0 || x == 1;
addOptional(p,'flag',0,fcn);
addParameter(p,'var',"r-");   
addParameter(p,'fillvar',"r");      % 默认参数
parse(p,varargin{:});
x1 = c + s * [cos(th) sin(th);-sin(th) cos(th)];
x2 = c + [-s(1) s(2)] * [cos(th) sin(th);-sin(th) cos(th)];
x3 = c + [-s(1) -s(2)] * [cos(th) sin(th);-sin(th) cos(th)];
x4 = c + [s(1) -s(2)] * [cos(th) sin(th);-sin(th) cos(th)];
plot([x1(1) x2(1) x3(1) x4(1) x1(1)], [x1(2) x2(2) x3(2) x4(2) x1(2)],p.Results.var);
if(p.Results.flag)
    fill([x1(1) x2(1) x3(1) x4(1) x1(1)], [x1(2) x2(2) x3(2) x4(2) x1(2)],p.Results.fillvar);
end
end
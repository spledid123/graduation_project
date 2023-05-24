function [] = plot_cir(x0,y0,r0,varargin)
%   [x0,y0]为圆心,r0为半径，画圆
%   flag为可选参数,flag==0不填充,flag==1填充
%   名称-值
%   var为字符串，plot的默认格式为'r-',为画圆的参数
%   linewidth为plot的线条宽度，默认2磅
%   fillvar为fill的颜色，默认'r';
%   示例：plot_cir(0,0,100,1,'var','r--','fillvar','r');plot_cir(0,0,100);
p = inputParser;            % 函数的输入解析器
fcn = @(x) x == 0 || x == 1;
addOptional(p,'flag',0,fcn);
addParameter(p,'var',"r-");      % 默认参数
addParameter(p,'fillvar',"r");      % 默认参数
addParameter(p,'linewidth',2);      % 默认参数
parse(p,varargin{:});       % 对输入变量进行解析，如果检测到前面的变量被赋值，则更新变量取值
th=0:0.1:2*pi+0.1;
x=x0+r0*cos(th);
y=y0+r0*sin(th);
plot(x,y,p.Results.var,'LineWidth',p.Results.linewidth);
if(p.Results.flag)
    fill(x,y,p.Results.fillvar,'LineStyle','none');
end
end
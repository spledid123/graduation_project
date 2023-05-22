%   针对六边形孔喉结构
%   计算在孔/喉内弹道长度的分布
%   计算孔/喉单次碰撞离开的概率
%   
clc;
clear;
%   输入参数
R = 1.1;% （喉宽度+孔直径）/孔直径
a = 0.25;% 喉道长度/孔直径

wid = 2*(R-1);% 喉道宽度
t2 = (4*(R-1)^2+4*a^2)^0.5;%    有限喉道的最长弹道长度
fv = readtable('data\pdf\Bol_2D_pdf_T_40.txt');%    速度分布表[v,f(v)]
filename1 = 'data\pdf\pdf_throat_len_R_';
filename_R = num2str(R);
filename_a = num2str(a);
filename2 = '_a_';
filename3 = '.txt';
filename = strcat(filename1,filename_R,filename2,filename_a,filename3);
fp = readtable(filename);%    喉道弹道长度概率分布表[r,f(r)]
fp = table2array(fp);
fp = fp(:,2);
fv = table2array(fv);
fv = fv(:,2);
%%   孔的弹道长度分布
ls = 2e-5/1e2;
po_x = 0:ls:2e-5;
po_y=[];
for i=1:length(po_x)
    po_y(i) = po_pdf(po_x(i));
end
figure;
plot(po_x(1:end-1),po_y(1:end-1),'DisplayName','解析结果');
title('孔内弹道长度分布','fontsize',14);
xlabel('孔内弹道长度','fontsize',14);
ylabel('概率密度','FontSize',14);
%%   喉的弹道长度分布
ls = (t2-wid)/1e2;
th_x = wid*(1.001):ls:t2;
th_x = th_x ./ 1e5;
th_y=[];
th_yy=[];
for i=1:length(th_x)
    th_y(i) = th_len_pdf(th_x(i),R,a);
end
figure;
plot(th_x,th_y,'DisplayName','解析结果');
title('有限长度喉道弹道长度分布','fontsize',14);
%%   孔内每次反射的出射概率
Pp = 3*asin(R-1)/pi;% 用面积比例估计
%%   有限喉道每次反射出射的平均概率
tx = 1e4;
ls = a/tx;
x = 0:ls:a;
y=[];
for i=1:length(x)
    y(i)=fx4(x(i),R,a);
end
Pt = trapz(x,y);
%%   f(t)=f(r/v);
%   孔内的单位时间分布
s=[];
x = 0:1e-9:1e-5;
mu = 0;
for i=1:length(x)
%     s(i)=ft_po(x(i),40);
    s(i)=ft_po_exl(x(i),40,fv);
end
figure;
plot(x,s,'DisplayName','孔内单位时间分布');
Tp=trapz(x,s.*x);
xlabel('孔内单步时间','FontSize',14);
ylabel('概率密度','FontSize',14);
%%   事先做好值表，方便求积分
%   喉内单位时间分布
s=[];
x = 0:5e-10:1e-6;
mu = 0;
for i=1:length(x)
%     s(i) = ft_th(x(i),R,a,fv);
    s(i) = ft_th_exl(x(i),R,a,fp,fv);
end
figure;
plot(x,s);
Tt=trapz(x,s.*x);
%%  单格点时间
T0 = (1/Pp)*Tp+3/2*(1/Pt)*Tt;
f = f_2D(R,a,T0);







function [sd] = ft_th_exl(t,R,a,fp,fv)
%    求喉道弹道时间分布
    x1 = 0:1e-1:100;
    x2 = 100:1e-2:400;
    x3 = 400:1e-1:600;
    x=[x1 x2 x3];
    y=[];
    for i=1:length(x)
        y(i)=x(i)*th_len_pdf_exl(t*x(i),fp,R,a)*v_pdf_exl(x(i),fv);
    end
    sd = trapz(x,y);
end
function [sd] = ft_th(t,R,a,fv)
%    求喉道弹道时间分布
    x1 = 0:1e-1:100;
    x2 = 100:1e-2:400;
    x3 = 400:1e-1:600;
    x=[x1 x2 x3];
    for i=1:length(x)
        y(i)=x(i)*th_len_pdf(t*x(i),R,a)*v_pdf(x(i),40);
    end
    sd = trapz(x,y);
end
function sd = ft_po_exl(t,T,fV)
%    求孔的弹道时间分布
    x1 = 0:1e-1:100;
    x2 = 100:1e-2:400;
    x3 = 400:1e-1:600;
    x=[x1 x2 x3];
    y=[];
    for i=1:length(x)
        y(i)=x(i)*po_pdf(t*x(i))*v_pdf_exl(x(i),fV);
    end
    sd = trapz(x,y);
end
function sd = ft_po(t,T)
%    求孔的弹道时间分布
    x1 = 0:1e-1:100;
    x2 = 100:1e-2:400;
    x3 = 400:1e-1:600;
    x=[x1 x2 x3];
    y=[];
    for i=1:length(x)
        y(i)=x(i)*po_pdf(t*x(i))*v_pdf(x(i),T);
    end
    sd = trapz(x,y);
end
function sd = v_pdf_exl(v,fV)%  查表求f(v),T=40K
    k = floor(v/0.1)+1;
    if((k<length(fV))==0)
        sd = 0;
    else
        sd = (fV(k+1)-fV(k))/0.1*(v-0.1*(k-1))+fV(k);
    end
end
function sd = v_pdf(v,T)%   玻尔兹曼分布的pdf
    k = 1.380649e-23;
    NA = 6.02e23;
    m = 18e-3/NA;
    A = m/2/k/T;
    sd = 2*A*v*exp(-A*v^2);
end

function sd = po_pdf(r)%    圆形孔弹道长度分布
    r = r*1e5;
    if((r<2)==0)
        sd = 0;
    else
    sd = r/4*(1-r^2/4)^(-0.5);
    end
    sd = sd*1e5;
end
function sd = th_len_pdf_exl(v,fV,R,a)%   有限长喉道，f(r),查表插值,R=1.1,a=0.25
    v = v * 1e5;
    t1 = 2*(R-1);
    t2 = (4*(R-1)^2+4*a^2)^0.5;
    k = floor((v-t1)/0.0001)+1;
    if(v<t1)
        sd=0;
        return;
    end
    if(v>t2)
        sd=0;
    end
    if((k<length(fV))==0)
        sd = 0;
    else
    sd = (fV(k+1)-fV(k))/0.0001*(v-(t1+0.0001*(k-1)))+fV(k);
    end
    sd = sd*1e5;
end
function sd = th_len_pdf(r,R,a)%    有限长度的喉道弹道长度
%   2*(R-1)为喉道宽度，2*a为喉道长度，r为弹道长度
    r = r * 1e5;
    wid = 2*(R-1);% 喉道宽度；
    t2 = (wid^2+4*a^2)^0.5;
    c = 1/2/a*(t2-wid);
    if((r>wid)==0)
        sd=0;
    elseif(r>t2)
        sd=0;
    else
    sd = wid^2/r^2*(1/(r^2-wid^2)^0.5-1/2/a)/c;
    end
    sd = sd * 1e5;
end
function sd = fx4(x,R,a)
    sd = 1/2/a*(2-x/(4*(R-1)^2+x^2)^0.5-(2*a-x)/(4*(R-1)^2+(2*a-x)^2)^0.5);
end
function sd = th_pdf(r,R)%  无限长度喉道弹道长度
    sd = 4*(R-1)^2/r^3*(1-4*(R-1)^2/r^2)^(-0.5);
end
function sd = f_2D(R,a,meanDT)
r = 1;
t = (1 + a) * r;%   六边形边长的一半
h = 2*r*(R-1);
len = 2*a*r;
A = pi + 1.5 * h * len;%    单孔面积
A = A * 1e-10;%    1m:10微米
h = h * 1e-5;
k = 1.380649e-23;
NA = 6.02e23;
m = 18e-3/NA;
v = (pi*k*40/2/m)^0.5;
sd = 3*v*meanDT*h/A/pi;
end
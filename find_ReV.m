%   圆阵列
%   给定一个点，以它为圆心画圆，计算孔隙度随半径的变化。得到孔隙结构的ReV
R = 1.2;
Rmax = 1;
Rmin = 1;
kR = [];
x0 = 0;
y0 = 0;
%   随机生成圆阵列
for i = -100:100
    for j = -100:100
        kR(fun_fromktoflag(i,j)) = rand()*(Rmax - Rmin)+Rmin;
    end
end
phi = [];
r0 = 0:0.12:50;
for i = 1:length(r0)
    S0 = pi*r0(i)^2;
    S = 0;% 固体面积
    for m = -100:100
    for n = -100:100
        S = S + fun_get_s(x0,y0,r0(i),(2*m+1)*R,(2*n+1)*R,kR(fun_fromktoflag(m,n)));
    end
    end
    phi(i) = 1-S/S0;
end
plot(pi*r0.^2,phi);
%%  找到粒子数密度的ReV
clc;
clear;
R = 1.2;
Rmax = 1;
Rmin = 1;
filename = 'data\onepartical_long_path_cir\rxT_circle_onepartical_cubecir_R_1.2_rmax_1.1_rmin_0.9(1).txt';
A = readtable(filename);
X = table2array(A(:,2));
Y = table2array(A(:,3));
tX = [];
tY = [];
kx = [];
ky = [];
flag  = [];
sd = [];
meansd = [];
si = [];
rx=[];
TH = [];
hold on;
% set(gca,'xscale','log');
for t = 1:100
    sd = [];
    for i = 1:length(X)
        [tX(i),tY(i),kx(i),ky(i),rx(i)] = fun_toPBC(X(i),Y(i),R,t);
    end
    j = 1;
    for i = 2:length(X)
        if(kx(i)~=kx(i-1) || ky(i)~=ky(i-1))
            flag(j) = i - 1;
            j = j + 1;
        end
    end
    flag = flag(1:j-1);
    meansd = 0;
    j=1;
    for i = 1:length(flag)-1
        sd(i) = fun_get_len(flag(i),flag(i+1),X,Y);
    end
    for i = 1:length(flag)-1
        if(kx(flag(i))~=kx(flag(i)+1) && ky(flag(i))==ky(flag(i)+1))
        TH(j) = fun_get_TH(flag(i),flag(i) + 1,X,Y);
        j=j+1;
        end
    end
    for i = 1:length(sd)
        meansd(i) = sum(sd(1:i))/i;
    end
    n = 10;
kA = R^2-pi*((Rmax - Rmin)^2/12-((Rmax+Rmin)/2)^2);
    i = 1:length(sd);
    ReV = (i*kA.^0.5./n).^1.5;
%     plot(ReV,meansd);
si(t) = sum(sd)/length(sd);
% %     plot(t,sum(meansd),'o');
end
t = 1:20;
plot(t,si);
% p = polyfit(t,si,1);
% p = p(1);
% n = 10;
% kA = R^2-pi*((Rmax - Rmin)^2/12-((Rmax+Rmin)/2)^2);















function x = fun_fromktoflag(kx,ky)
    %得到kR的序号
    % (kx,ky)指圆心在((2*kx+1)*R,(2*ky+1)*R)的圆。
    if(kx==0 && ky==0)
        x=1;
        return;
    end
    t = abs(kx)+abs(ky);
    x = 2*t*(t-1)+1;
    if(kx == 0)
        if(ky > 0)
            x = x + ky + 1;
        else
            x = x - ky * 3 + 1;
        end
    elseif(ky == 0)
        if(kx > 0)
            x = x + 1;
        else
            x = x - kx * 2 + 1;
        end
    else
        if(kx > 0 && ky > 0)
            x = x + ky + 1;
        elseif(kx < 0 && ky > 0)
            x = x + t - kx + 1;
        elseif(kx < 0 && ky < 0)
            x = x + t * 2 - ky + 1;
        else
            x = x + t * 3 + kx + 1;
        end
    end
end
function sd = fun_get_s(x0,y0,r0,x,y,r)
%   以(x0,y0)画半径为r0的圆，(x,y)为圆心，半径为r的圆在其中的面积为sd
t = ((x0-x)^2+(y0-y)^2)^0.5;
if(t > (r0 + r))
    sd = 0;
    return;
elseif(t > abs(r0 - r))
    %   扇形减去三角形的面积
    alp = acos((r0^2 + t^2 -r^2)/2/r0/t);
    alp1 = acos((r^2 + t^2 -r0^2)/2/r/t);
    sd = alp * r0^2 - r0^2*sin(alp)*cos(alp) + alp1 * r^2 - r^2*sin(alp1)*cos(alp1);
else
    if(r0 < r)
        sd = pi*r0^2;
    else
        sd = pi*r^2;
    end
end
end
function [x,y,kx,ky,rx] = fun_toPBC(x0,y0,R,t)
%   将非周期性边界得到的结果化为周期性边界
%   以四个圆心围成的方格为一格，以(0,0)所在格点为中心，t表示周期的小方块的边长，如果它是计数，(0,0)在中心；如果是偶数
%   （0，0）小格在中心线偏左偏上
a = 2 * R;% 小格的边长
if(mod(t,2) == 1)
    kx = floor((x0 - t*a/2)/t/a);
    ky = floor((y0 - t*a/2)/t/a);
    x = x0 - t*a*kx - t*a/2;
    y = y0 - t*a*ky - t*a/2;
else
    kx = floor((x0 - (t+1)*a/2)/t/a);
    ky = floor((y0 - (t+1)*a/2)/t/a);
    x = x0 - t*a*kx - (t+1)*a/2;
    y = y0 - t*a*ky - (t+1)*a/2;
end
if(mod(kx,3)==0)
    if(mod(ky,3) == 0)
        rx = 5;
    elseif(mod(ky,3)==1)
        rx = 2;
    else
        rx = 8;
    end
elseif(mod(kx,3)==1)
    if(mod(ky,3) == 0)
        rx = 6;
    elseif(mod(ky,3)==1)
        rx = 3;
    else
        rx = 9;
    end
else
        if(mod(ky,3) == 0)
        rx = 4;
    elseif(mod(ky,3)==1)
        rx = 1;
    else
        rx = 7;
    end
end
end
function sd = fun_get_len(i,j,X,Y)
if((j > i) == 0)
    sd = 0;
    error('ji');
end
sd = 0;
for t = i + 1 : j
    sd = sd + ((X(t)-X(t-1))^2+(Y(t)-Y(t-1))^2)^0.5;
end
end
function th = fun_get_TH(i,j,X,Y)
if((j > i) == 0)
    sd = 0;
    error('ji');
end
th = atan2(Y(j)-Y(i),X(j)-X(i));
% if(abs(th)>pi/4 && abs(th)<3*pi/4)
%     th = th - pi/2;
% end
if(th>pi/2)
    th = th-pi;
elseif(th<-pi/2)
    th=th+pi;
end
end
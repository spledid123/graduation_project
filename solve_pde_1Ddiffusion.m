% 解一维的扩散偏微分方程，有限元

%   边界条件：无限长
% c0=0.0027;
% t=10*60;
% D=1.52;%cm
% fun=@(x)c0*erfc(x/2/(D*t)^0.5);
% y=[];
% x=30:0.1:1000;
% len=length(x);
% for i=1:len
%     y(i)=fun(x(i));
% end
% s=[];
% for i=2:len
% s(i)=trapz(x(1:i),y(1:i));
% end
% tt=s(end)*pi*1^2/10000/4;

%   边界条件：有界
clc;
clear;
% 解网格
L = 12.8e-2;
x2 = 0.01;
tmax = 0.001;
tmaxn = 500;
xmaxn = 100;
x = [linspace(0,L,xmaxn)];
t = [linspace(0,tmax,tmaxn)];

m = 0;
sol = pdepe(m,@heatpde,@heatic,@bcfun,x,t);
u1 = sol(:,:,1);

% surf(x,t,u1);
% title('u_1(x,t)')
% xlabel('Distance x')
% ylabel('Time t')
uend = u1(:,50);
uend1 = u1(:,49);
du = (uend1-uend)./(L/xmaxn);
figure;
ux=u1(tmaxn,:);
plot(t,uend/1e-4);
T = 300;
v=(8*8.314*T/pi/0.018)^0.5;
D = v*2*0.85*1e-2/3;
for i=1:tmaxn
    sumu(i) = D*(tmax/tmaxn)*sum(du(1:i))/1e-4;
end
hold on;
figure;
plot(t,sumu);
figure;
hold on;
ut = u1(1,:);
utt = u1(end,:);
plot(x,ut);
plot(x,utt);


%   偏微分方程
function [c,f,s] = heatpde(x,t,u,dudx)
T = 300;
v=(8*8.314*T/pi/0.018)^0.5;
D = v*2*0.85*1e-2/3;
c = 1;
f = D*dudx;
s = 0;
end
%   初始条件
function u0 = heatic(x)
u0 = 0;
end
%   边界条件
function [pL,qL,pR,qR] = bcfun(xL,uL,xR,uR,t)
c0 = 1e-4;
pR = uR;
qR = 0;
% pR = 0;
% qR = 1;
pL = uL - c0;
qL = 0;
end



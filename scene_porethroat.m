function sd  = scene_porethroat(x, R, a)
%   sd(1):返回点到阵列的最短距离,sd(2),sd(3)为点最靠近的圆心的阵列编号(x,y方向);sd(4)表示其是否为孔，是为1，否为0
%   x=[x,y]为坐标
% R为喉道长度/直径；a为喉道宽度/直径

r = 1;  %   圆半径
t = (1 + a) * r;%   六边形边长的一半
l = 1.7320508075688772935274463415059;% 3^0.5
ky = floor(x(2)/l/t);
rymin = ky * l * t;
rymax = (ky + 1) * l * t;
kx = floor(x(1)/t);
kkx = floor(kx / 6);
%   由于按圆心分块，所以如果喉道过宽，跨过kx的分界线，就需要修正
if(mod(ky,2) == 1)
    if(mod(kx, 6)==1)
        rxmax = (kkx * 6 + 1) * t;
        if(atan2(x(2)-rymin,x(1)-rxmax)>pi/3)
            kx = kx - 1;
        end
    elseif(mod(kx, 6)==2)
        rxmin = (kkx * 6 + 3) * t;
        if(atan2(x(2)-rymin,x(1)-rxmin)<2*pi/3)
            kx = kx + 1;
        end
    end
        if(mod(kx, 6)== 4)
        rxmin = (kkx * 6 + 4) * t;
        if(atan2(x(2)-rymax,x(1)-rxmin)<-pi/3)
            kx = kx - 1;
        end
    elseif(mod(kx, 6)==5)
        rxmax = (kkx * 6 + 6) * t;
        if(atan2(x(2)-rymax,x(1)-rxmax)>-2*pi/3)
            kx = kx + 1;
        end
    end
elseif(mod(ky,2) == 0)
    if(mod(kx, 6)== 4)
        rxmin = (kkx * 6 + 4) * t;
        if(atan2(x(2)-rymin,x(1)-rxmin)>pi/3)
            kx = kx - 1;
        end
    elseif(mod(kx, 6)==5)
        rxmax = (kkx * 6 + 6) * t;
        if(atan2(x(2)-rymin,x(1)-rxmax)<2*pi/3)
            kx = kx + 1;
        end
    end
        if(mod(kx, 6)==1)
        rxmax = (kkx * 6 + 1) * t;
        if(atan2(x(2)-rymax,x(1)-rxmax)<-pi/3)
            kx = kx - 1;
        end
    elseif(mod(kx, 6)==2)
        rxmin = (kkx * 6 + 3) * t;
        if(atan2(x(2)-rymax,x(1)-rxmin)>-2*pi/3)
            kx = kx + 1;
        end
    end
end
kkx = floor(kx / 6);
%   如此分块就可以避免计算其它两个喉道的距离，于是只考虑两个圆以及喉道，这就是twocircleSDF的原理
sd = [1 1 1 1];
if(mod(ky,2)==0)%   ky是偶数
    if(mod(kx,6)==4||mod(kx,6)==5)
        rxmin = (kkx * 6 + 4) * t;
        rxmax = (kkx * 6 + 6) * t;
        [sd(1),sd(2),sd(3),sd(4)] = twocircleSDF(x, rxmin,rymin,rxmax,rymin,r,2*r*(R(1)-1), 2*a*r, 0);
    elseif(mod(kx,6) == 0)
        rxmin = (kkx * 6) * t;
        rxmax = (kkx * 6 + 1) * t;
        [sd(1),sd(2),sd(3),sd(4)] = twocircleSDF(x, rxmin,rymin,rxmax,rymax,r,2*r*(R(1)-1), 2*a*r, 1);
    elseif(mod(kx,6)==1||mod(kx,6)==2)
        rxmin = (kkx * 6 + 1) * t;
        rxmax = (kkx * 6 + 3) * t;  
        [sd(1),sd(2),sd(3),sd(4)] = twocircleSDF(x, rxmin,rymax,rxmax,rymax,r,2*r*(R(1)-1), 2*a*r, 0);
    else
        rxmin = (kkx * 6 + 3) * t;
        rxmax = (kkx * 6 + 4) * t;
        [sd(1),sd(2),sd(3),sd(4)] = twocircleSDF(x, rxmin,rymax,rxmax,rymin,r,2*r*(R(1)-1), 2*a*r, -1);
    end
else
    if(mod(kx,6)==1||mod(kx,6)==2)
        rxmin = (kkx * 6 + 1) * t;
        rxmax = (kkx * 6 + 3) * t;
        [sd(1),sd(2),sd(3),sd(4)] = twocircleSDF(x, rxmin,rymin,rxmax,rymin,r,2*r*(R(1)-1), 2*a*r, 0);
    elseif(mod(kx,6) == 3)
        rxmin = (kkx * 6 + 3) * t;
        rxmax = (kkx * 6 + 4) * t;
        [sd(1),sd(2),sd(3),sd(4)] = twocircleSDF(x, rxmin,rymin,rxmax,rymax,r,2*r*(R(1)-1), 2*a*r, 1);
    elseif(mod(kx,6)==4||mod(kx,6)==5)
        rxmin = (kkx * 6 + 4) * t;
        rxmax = (kkx * 6 + 6) * t;
        [sd(1),sd(2),sd(3),sd(4)] = twocircleSDF(x, rxmin,rymax,rxmax,rymax,r,2*r*(R(1)-1), 2*a*r, 0);
    else        
        rxmin = (kkx * 6) * t;
        rxmax = (kkx * 6 + 1) * t;
        [sd(1),sd(2),sd(3),sd(4)] = twocircleSDF(x, rxmin,rymax,rxmax,rymin,r,2*r*(R(1)-1), 2*a*r, -1);
    end
end
end
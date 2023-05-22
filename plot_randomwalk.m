%   画示意图——六边形孔喉随机游走模型
hold on;
axis equal;
axis off;
set(gca,'xlim',[-3.5 3.3]);
plot_cir(0,0,1);
plot_random_2(cos(pi/3),sin(pi/3),pi/3);
plot_random_2(cos(-pi),sin(-pi),-pi);
plot_random_2(cos(-pi/3),sin(-pi/3),-pi/3);
text(0,0,'孔','FontSize',14);
text(-2.2,0.5,'喉','FontSize',14);
title('(a) 2元随机游走模型','FontSize',14);



function [] = plot_random_2(x,y,theta)
x1 = x+0.7*cos(theta);
y1 = y+0.7*sin(theta);
plot([x x1],[y y1],'r');
th = theta-pi/2;
x2 = x1+0.25*cos(th);
y2 = y1+0.25*sin(th);
x3 = x1-0.25*cos(th);
y3 = y1-0.25*sin(th);
plot([x1 x2],[y1 y2],'r');
plot([x1 x3],[y1 y3],'r');
x4 = x2+1*cos(theta);
y4 = y2+1*sin(theta);
plot([x2 x4],[y2 y4],'r');
x5 = x3+1*cos(theta);
y5 = y3+1*sin(theta);
plot([x3 x5],[y3 y5],'r');
plot([x4 x5],[y4 y5],'r');
x6 = x1+1*cos(theta);
y6 = y1+1*sin(theta);
x7 = x6+0.7*cos(theta);
y7 = y6+0.7*sin(theta);
plot([x6 x7],[y6 y7],'r');
plot([(x2+x4)/2+0.2*cos(th) (x3+x5)/2-0.2*cos(th)],[(y2+y4)/2+0.2*sin(th) (y3+y5)/2-0.2*sin(th)],'r--');
end
function sd = scene(x, R, varargin)
%   返回点到边界的最短距离
%   x=[x,y]为坐标
%   R=[Rx,Ry]为多孔介质形状参数
p = inputParser;            % 函数的输入解析器
addOptional(p,'pm',0);
parse(p,varargin{:});
%%   圆形阵列
%   圆的半径为r,于是x方向孔的半径为Rx-1;y方向孔的半径为Ry-1;
%   于是x方向圆心坐标为(Rx)*(2*kx+1);y方向圆心坐标为Ry*(2*ky+1);
%   R = [Rx,Ry,半径r]
% if(length(R) ~= 3)
%     error('多孔介质的参数错了,scene.m错了');
% end
% R(1) = R(1) * R(3);
% R(2) = R(2) * R(3);
% kx = floor((x(1)/R(1)-1)/2);
% ky = floor((x(2)/R(2)-1)/2);
% 
% rxmin = R(1) * (2 * kx + 1);
% rxmax = R(1) * (2 * kx + 3);
% rymin = R(2) * (2 * ky + 1);
% rymax = R(2) * (2 * ky + 3);
% 
% a=[];
% a(1) = circleSDF(x,[rxmin,rymin],R(3));
% a(2) = circleSDF(x,[rxmax,rymin],R(3));
% a(3) = circleSDF(x,[rxmin,rymax],R(3));
% a(4) = circleSDF(x,[rxmax,rymax],R(3));
% 
% sd = min(a);

%%   孔喉六边形结构
% R(1)为喉道长度/直径；R(2)为喉道宽度/直径
% r = 1;  %   圆半径
% t = (1 + R(2)) * r;%   六边形边长的一半
% l = 1.7320508075688772935274463415059;% 3^0.5
% ky = floor(x(2)/l/t);
% rymin = ky * l * t;
% rymax = (ky + 1) * l * t;
% kx = floor(x(1)/t);
% kkx = floor(kx / 6);
% %   由于按圆心分块，所以如果喉道过宽，跨过kx的分界线，就需要修正
% if(mod(ky,2) == 1)
%     if(mod(kx, 6)==1)
%         rxmax = (kkx * 6 + 1) * t;
%         if(atan2(x(2)-rymin,x(1)-rxmax)>pi/3)
%             kx = kx - 1;
%         end
%     elseif(mod(kx, 6)==2)
%         rxmin = (kkx * 6 + 3) * t;
%         if(atan2(x(2)-rymin,x(1)-rxmin)<2*pi/3)
%             kx = kx + 1;
%         end
%     end
%         if(mod(kx, 6)== 4)
%         rxmin = (kkx * 6 + 4) * t;
%         if(atan2(x(2)-rymax,x(1)-rxmin)<-pi/3)
%             kx = kx - 1;
%         end
%     elseif(mod(kx, 6)==5)
%         rxmax = (kkx * 6 + 6) * t;
%         if(atan2(x(2)-rymax,x(1)-rxmax)>-2*pi/3)
%             kx = kx + 1;
%         end
%     end
% elseif(mod(ky,2) == 0)
%     if(mod(kx, 6)== 4)
%         rxmin = (kkx * 6 + 4) * t;
%         if(atan2(x(2)-rymin,x(1)-rxmin)>pi/3)
%             kx = kx - 1;
%         end
%     elseif(mod(kx, 6)==5)
%         rxmax = (kkx * 6 + 6) * t;
%         if(atan2(x(2)-rymin,x(1)-rxmax)<2*pi/3)
%             kx = kx + 1;
%         end
%     end
%         if(mod(kx, 6)==1)
%         rxmax = (kkx * 6 + 1) * t;
%         if(atan2(x(2)-rymax,x(1)-rxmax)<-pi/3)
%             kx = kx - 1;
%         end
%     elseif(mod(kx, 6)==2)
%         rxmin = (kkx * 6 + 3) * t;
%         if(atan2(x(2)-rymax,x(1)-rxmin)>-2*pi/3)
%             kx = kx + 1;
%         end
%     end
% end
% kkx = floor(kx / 6);
% %   如此分块就可以避免计算其它两个喉道的距离，于是只考虑两个圆以及喉道，这就是twocircleSDF的原理
% sd = [1 1 1 1];
% if(mod(ky,2)==0)%   ky是偶数
%     if(mod(kx,6)==4||mod(kx,6)==5)
%         rxmin = (kkx * 6 + 4) * t;
%         rxmax = (kkx * 6 + 6) * t;
%         [sd(1),sd(2),sd(3),sd(4)] = twocircleSDF(x, rxmin,rymin,rxmax,rymin,r,2*r*(R(1)-1), 2*R(2)*r, 0);
%     elseif(mod(kx,6) == 0)
%         rxmin = (kkx * 6) * t;
%         rxmax = (kkx * 6 + 1) * t;
%         [sd(1),sd(2),sd(3),sd(4)] = twocircleSDF(x, rxmin,rymin,rxmax,rymax,r,2*r*(R(1)-1), 2*R(2)*r, 1);
%     elseif(mod(kx,6)==1||mod(kx,6)==2)
%         rxmin = (kkx * 6 + 1) * t;
%         rxmax = (kkx * 6 + 3) * t;
%         [sd(1),sd(2),sd(3),sd(4)] = twocircleSDF(x, rxmin,rymax,rxmax,rymax,r,2*r*(R(1)-1), 2*R(2)*r, 0);
%     else
%         rxmin = (kkx * 6 + 3) * t;
%         rxmax = (kkx * 6 + 4) * t;
%         [sd(1),sd(2),sd(3),sd(4)] = twocircleSDF(x, rxmin,rymax,rxmax,rymin,r,2*r*(R(1)-1), 2*R(2)*r, -1);
%     end
% else
%     if(mod(kx,6)==1||mod(kx,6)==2)
%         rxmin = (kkx * 6 + 1) * t;
%         rxmax = (kkx * 6 + 3) * t;
%         [sd(1),sd(2),sd(3),sd(4)] = twocircleSDF(x, rxmin,rymin,rxmax,rymin,r,2*r*(R(1)-1), 2*R(2)*r, 0);
%     elseif(mod(kx,6) == 3)
%         rxmin = (kkx * 6 + 3) * t;
%         rxmax = (kkx * 6 + 4) * t;
%         [sd(1),sd(2),sd(3),sd(4)] = twocircleSDF(x, rxmin,rymin,rxmax,rymax,r,2*r*(R(1)-1), 2*R(2)*r, 1);
%     elseif(mod(kx,6)==4||mod(kx,6)==5)
%         rxmin = (kkx * 6 + 4) * t;
%         rxmax = (kkx * 6 + 6) * t;
%         [sd(1),sd(2),sd(3),sd(4)] = twocircleSDF(x, rxmin,rymax,rxmax,rymax,r,2*r*(R(1)-1), 2*R(2)*r, 0);
%     else
%         rxmin = (kkx * 6) * t;
%         rxmax = (kkx * 6 + 1) * t;
%         [sd(1),sd(2),sd(3),sd(4)] = twocircleSDF(x, rxmin,rymax,rxmax,rymin,r,2*r*(R(1)-1), 2*R(2)*r, -1);
%     end
% end
% sd = sd(1);
%%  周期性空腔圆阵列交错
%   整个周期性边界大致为九宫格：存在两种相同的矩形区块，分别有不同的形态，任意一种的上下左右都是另一种
%   圆阵列为numx*numy
% numx = 20;
% numy = 20;
% tx = 2 * numx * R(1);
% ty = 2 * numy * R(2);
% kx = floor(x(1)/tx);
% ky = floor(x(2)/ty);
% %   归一
% if(mod(kx,2) == mod(ky,2))  %   圆阵列
%     x(1) = x(1) - kx * tx;
%     x(2) = x(2) - ky * ty;
%     kx = floor((x(1)/R(1)-1)/2);%   在圆阵列的序号，左下角圆心为(0,0)
%     ky = floor((x(2)/R(2)-1)/2);
%     rxmin = R(1) * (2 * kx + 1);
%     rxmax = R(1) * (2 * kx + 3);
%     rymin = R(2) * (2 * ky + 1);
%     rymax = R(2) * (2 * ky + 3);
%     a=[];
%     a(1) = circleSDF(x,[rxmin,rymin],1);
%     a(2) = circleSDF(x,[rxmax,rymin],1);
%     a(3) = circleSDF(x,[rxmin,rymax],1);
%     a(4) = circleSDF(x,[rxmax,rymax],1);
%     %   边缘的情况
%     t = max(tx,ty);%    大值
%     if(kx < 0)% rxmin的圆不能要了（不存在）
%         a(1) = t;
%         a(3) = t;
%     elseif(kx > (numx - 2))% rxmax的圆不能要了（不存在）
%         a(2) = t;
%         a(4) = t;
%     end
%     if(ky < 0)% rymin的圆不能要了（不存在）
%         a(1) = t;
%         a(2) = t;
%     elseif(ky > (numy - 2))% rymax的圆不能要了（不存在）
%         a(3) = t;
%         a(4) = t;
%     end
%     sd = min(a);
% else    %   空腔，四个顶点为(0,0),(tx,0),(0,ty),(tx,ty)
%     x(1) = x(1) - kx * tx;
%     x(2) = x(2) - ky * ty;
%     sd = min([x(1),tx-x(1),x(2),ty-x(2)]);
%     sd = sd + R(1) - 1;%    多算一点，更快跨越空腔边界
% end
%%   整个周期性边界大致为九宫格：存在两种边长相等正方形区块，分别有不同的形态，任意一种的上下左右都是另一种
%   为了对称性，先暂定两种区块都是均质圆阵列且xy方向性质一样
%   多孔介质参数R=[右边的圆半径 右边的圆半径+喉道宽度一半 右边的一排/列圆数量 左边的圆半径 左边的圆半径+喉道宽度一半 左边的一排/列圆数量]
% if(length(R) ~= 6)
%     error('多孔介质的参数错了,scene.m错了');
% end
% numr = R(3);    % 右边圆阵列圆数量
% numl = R(6);    % 左边圆阵列圆数量
% Rr = R(1);      % 右边圆阵列圆半径
% Rl = R(4);      % 左边圆阵列圆半径
% tr = R(2);      % 右边圆阵列圆半径+喉道宽度一半
% tl = R(5);      % 左边圆阵列圆半径+喉道宽度一半
% if((tr < Rr) || (tl < Rl))
%     error("圆阵列的圆挤在一块了,多孔介质参数设置有问题");
% end
% if(abs(numr * tr - numl * tl) > 1e-5)
%     error("两个区块的边长不同,多孔介质参数设置有问题");
% end
% t = 2 * numr * tr;%  正方形区块边长
% kx = floor(x(1)/t);
% ky = floor(x(2)/t);
% x(1) = x(1) - kx * t;% 归一到左下端为(0,0)的区块
% x(2) = x(2) - ky * t;
% %   归一并求SDF
% if(mod(kx,2) == mod(ky,2))  %   右边圆阵列
%     kx = floor((x(1)/tr-1)/2);%   在圆阵列的序号，左下角圆心为(0,0)
%     ky = floor((x(2)/tr-1)/2);
%     rxmin = tr * (2 * kx + 1);
%     rxmax = tr * (2 * kx + 3);
%     rymin = tr * (2 * ky + 1);
%     rymax = tr * (2 * ky + 3);
%     a=[];
%     a(1) = circleSDF(x,[rxmin,rymin],Rr);
%     a(2) = circleSDF(x,[rxmax,rymin],Rr);
%     a(3) = circleSDF(x,[rxmin,rymax],Rr);
%     a(4) = circleSDF(x,[rxmax,rymax],Rr);
%     a(5) = t;%  可能的临近其他区块的圆
%     a(6) = t;
%     a(7) = t;
%     a(8) = t;
%     %   边缘的情况
%     if(kx < 0)% rxmin的圆不能要了（不存在）
%         a(1) = t;
%         a(3) = t;
%         rxmin = -tl;
%         kky = floor((x(2)/tl-1)/2);
%         rymin = tl * (2 * kky + 1);
%         rymax = tl * (2 * kky + 3);
%         a(5) = circleSDF(x,[rxmin,rymin],Rl);
%         a(6) = circleSDF(x,[rxmin,rymax],Rl);
%         if(kky < 0)
%             a(5) = t;
%         elseif(kky > (numl - 2))
%             a(6) = t;
%         end
%     elseif(kx > (numr - 2))% rxmax的圆不能要了（不存在）
%         a(2) = t;
%         a(4) = t;
%         rxmax = t + tl;
%         kky = floor((x(2)/tl-1)/2);
%         rymin = tl * (2 * kky + 1);
%         rymax = tl * (2 * kky + 3);
%         a(5) = circleSDF(x,[rxmax,rymin],Rl);
%         a(6) = circleSDF(x,[rxmax,rymax],Rl);
%         if(kky < 0)
%             a(5) = t;
%         elseif(kky > (numl - 2))
%             a(6) = t;
%         end
%     end
%     if(ky < 0)% rymin的圆不能要了（不存在）
%         a(1) = t;
%         a(2) = t;
%         rymin = -tl;
%         kkx = floor((x(1)/tl-1)/2);
%         rxmin = tl * (2 * kkx + 1);
%         rxmax = tl * (2 * kkx + 3);
%         a(7) = circleSDF(x,[rxmin,rymin],Rl);
%         a(8) = circleSDF(x,[rxmax,rymin],Rl);
%         if(kkx < 0)
%             a(7) = t;
%         elseif(kkx > (numl - 2))
%             a(8) = t;
%         end
%     elseif(ky > (numr - 2))% rymax的圆不能要了（不存在）
%         a(3) = t;
%         a(4) = t;
%         rymax = t + tl;
%         kkx = floor((x(1)/tl-1)/2);
%         rxmin = tl * (2 * kkx + 1);
%         rxmax = tl * (2 * kkx + 3);
%         a(7) = circleSDF(x,[rxmin,rymax],Rl);
%         a(8) = circleSDF(x,[rxmax,rymax],Rl);
%         if(kkx < 0)
%             a(7) = t;
%         elseif(kkx > (numl - 2))
%             a(8) = t;
%         end
%     end
%     sd = min(a);
% else    %   左边圆阵列
%     kx = floor((x(1)/tl-1)/2);%   在圆阵列的序号，左下角圆心为(0,0)
%     ky = floor((x(2)/tl-1)/2);
%     rxmin = tl * (2 * kx + 1);
%     rxmax = tl * (2 * kx + 3);
%     rymin = tl * (2 * ky + 1);
%     rymax = tl * (2 * ky + 3);
%     a=[];
%     a(1) = circleSDF(x,[rxmin,rymin],Rl);
%     a(2) = circleSDF(x,[rxmax,rymin],Rl);
%     a(3) = circleSDF(x,[rxmin,rymax],Rl);
%     a(4) = circleSDF(x,[rxmax,rymax],Rl);
%     a(5) = t;
%     a(6) = t;
%     a(7) = t;
%     a(8) = t;
%     %   边缘的情况
%     if(kx < 0)% rxmin的圆不能要了（不存在）
%         a(1) = t;
%         a(3) = t;
%         rxmin = -tr;
%         kky = floor((x(2)/tr-1)/2);
%         rymin = tr * (2 * kky + 1);
%         rymax = tr * (2 * kky + 3);
%         a(5) = circleSDF(x,[rxmin,rymin],Rr);
%         a(6) = circleSDF(x,[rxmin,rymax],Rr);
%         if(kky < 0)
%             a(5) = t;
%         elseif(kky > (numr - 2))
%             a(6) = t;
%         end
%     elseif(kx > (numl - 2))% rxmax的圆不能要了（不存在）
%         a(2) = t;
%         a(4) = t;
%         rxmax = t + tr;
%         kky = floor((x(2)/tr-1)/2);
%         rymin = tr * (2 * kky + 1);
%         rymax = tr * (2 * kky + 3);
%         a(5) = circleSDF(x,[rxmax,rymin],Rr);
%         a(6) = circleSDF(x,[rxmax,rymax],Rr);
%         if(kky < 0)
%             a(5) = t;
%         elseif(kky > (numr - 2))
%             a(6) = t;
%         end
%     end
%     if(ky < 0)% rymin的圆不能要了（不存在）
%         a(1) = t;
%         a(2) = t;
%         rymin = -tr;
%         kkx = floor((x(1)/tr-1)/2);
%         rxmin = tr * (2 * kkx + 1);
%         rxmax = tr * (2 * kkx + 3);
%         a(7) = circleSDF(x,[rxmin,rymin],Rr);
%         a(8) = circleSDF(x,[rxmax,rymin],Rr);
%         if(kkx < 0)
%             a(7) = t;
%         elseif(kkx > (numr - 2))
%             a(8) = t;
%         end
%     elseif(ky > (numl - 2))% rymax的圆不能要了（不存在）
%         a(3) = t;
%         a(4) = t;
%         rymax = t + tr;
%         kkx = floor((x(1)/tr-1)/2);
%         rxmin = tr * (2 * kkx + 1);
%         rxmax = tr * (2 * kkx + 3);
%         a(7) = circleSDF(x,[rxmin,rymax],Rr);
%         a(8) = circleSDF(x,[rxmax,rymax],Rr);
%         if(kkx < 0)
%             a(7) = t;
%         elseif(kkx > (numr - 2))
%             a(8) = t;
%         end
%     end
%     sd = min(a);
% end
%%   整个周期性边界大致为九宫格：存在两种边长相等正方形区块，分别有不同的形态，任意一种的上下左右都是另一种，且之间存在空余消除边界效应
%   为了对称性，先暂定两种区块都是均质圆阵列且xy方向性质一样
%   为了消除边界效应，四周空出一定距离
%   多孔介质参数R=[右边的圆半径 右边的圆半径+喉道宽度一半 右边的一排/列圆数量 左边的圆半径 左边的圆半径+喉道宽度一半 左边的一排/列圆数量 空出的长度]
% if(length(R) ~= 8)
%     error('多孔介质的参数错了,scene.m错了');
% end
% numr = R(3);    % 右边圆阵列圆数量
% numl = R(6);    % 左边圆阵列圆数量
% Rr = R(1);      % 右边圆阵列圆半径
% Rl = R(4);      % 左边圆阵列圆半径
% tr = R(2);      % 右边圆阵列圆半径+喉道宽度一半
% tl = R(5);      % 左边圆阵列圆半径+喉道宽度一半
% tsdr = R(7);     %  圆阵列与方格空出的距离
% tsdl = R(8);
% if((tr < Rr) || (tl < Rl))
%     error("圆阵列的圆挤在一块了,多孔介质参数设置有问题");
% end
% if(abs(numr * tr + tsdr - numl * tl - tsdl) > 1e-5)
%     error("两个区块的边长不同,多孔介质参数设置有问题");
% end
% t = 2 * numr * tr + tsdr * 2;%  正方形区块边长
% kx = floor(x(1)/t);
% ky = floor(x(2)/t);
% x(1) = x(1) - kx * t;% 归一到左下端为(0,0)的区块
% x(2) = x(2) - ky * t;
% %   归一并求SDF
% if(mod(kx,2) == mod(ky,2))  %   右边圆阵列
%     kx = floor(((x(1) - tsdr)/tr-1)/2);%   在圆阵列的序号，左下角圆心为(0,0)
%     ky = floor(((x(2) - tsdr)/tr-1)/2);
%     if(kx < 0) %    让计算出的圆心位置正确
%         kx = -1;
%     elseif(kx > (numr - 2))
%         kx = numr - 1;
%     end
%     if(ky < 0)
%         ky = -1;
%     elseif(ky > (numr - 2))
%         ky = numr - 1;
%     end
%     rxmin = tr * (2 * kx + 1) + tsdr;
%     rxmax = tr * (2 * kx + 3) + tsdr;
%     rymin = tr * (2 * ky + 1) + tsdr;
%     rymax = tr * (2 * ky + 3) + tsdr;
%     a=[];
%     a(1) = circleSDF(x,[rxmin,rymin],Rr);
%     a(2) = circleSDF(x,[rxmax,rymin],Rr);
%     a(3) = circleSDF(x,[rxmin,rymax],Rr);
%     a(4) = circleSDF(x,[rxmax,rymax],Rr);
%     a(5) = t;%  可能的临近其他区块的圆     
%     a(6) = t;
%     a(7) = t;
%     a(8) = t;
%     %   边缘的情况
%     if(kx < 0)% rxmin的圆不能要了（不存在）
%         a(1) = t;
%         a(3) = t;
%         rxmin = -tl - tsdl;
%         kky = floor(((x(2) - tsdl)/tl-1)/2);
%         rymin = tl * (2 * kky + 1) + tsdl;
%         rymax = tl * (2 * kky + 3) + tsdl;
%         a(5) = circleSDF(x,[rxmin,rymin],Rl);
%         a(6) = circleSDF(x,[rxmin,rymax],Rl);
%         if(kky < 0)
%             a(5) = t;
%         elseif(kky > (numl - 2))
%             a(6) = t;
%         end
%     elseif(kx > (numr - 2))% rxmax的圆不能要了（不存在）
%         a(2) = t;
%         a(4) = t;
%         rxmax = t + tl + tsdl;
%         kky = floor(((x(2) - tsdl)/tl-1)/2);
%         rymin = tl * (2 * kky + 1) + tsdl;
%         rymax = tl * (2 * kky + 3) + tsdl;
%         a(5) = circleSDF(x,[rxmax,rymin],Rl);
%         a(6) = circleSDF(x,[rxmax,rymax],Rl);
%         if(kky < 0)
%             a(5) = t;
%         elseif(kky > (numl - 2))
%             a(6) = t;
%         end
%     end
%     if(ky < 0)% rymin的圆不能要了（不存在）
%         a(1) = t;
%         a(2) = t;
%         rymin = -tl - tsdl;
%         kkx = floor(((x(1) - tsdl)/tl-1)/2);
%         rxmin = tl * (2 * kkx + 1) + tsdl;
%         rxmax = tl * (2 * kkx + 3) + tsdl;
%         a(7) = circleSDF(x,[rxmin,rymin],Rl);
%         a(8) = circleSDF(x,[rxmax,rymin],Rl);
%         if(kkx < 0)
%             a(7) = t;
%         elseif(kkx > (numl - 2))
%             a(8) = t;
%         end
%     elseif(ky > (numr - 2))% rymax的圆不能要了（不存在）
%         a(3) = t;
%         a(4) = t;
%         rymax = t + tl + tsdl;
%         kkx = floor(((x(1) - tsdl)/tl-1)/2);
%         rxmin = tl * (2 * kkx + 1) + tsdl;
%         rxmax = tl * (2 * kkx + 3) + tsdl;
%         a(7) = circleSDF(x,[rxmin,rymax],Rl);
%         a(8) = circleSDF(x,[rxmax,rymax],Rl);
%         if(kkx < 0)
%             a(7) = t;
%         elseif(kkx > (numl - 2))
%             a(8) = t;
%         end
%     end
%     sd = min(a);
% else    %   左边圆阵列
%     kx = floor(((x(1) - tsdl)/tl-1)/2);%   在圆阵列的序号，左下角圆心为(0,0)
%     ky = floor(((x(2) - tsdl)/tl-1)/2);
%     if(kx < 0) %    让计算出的圆心位置正确
%         kx = -1;
%     elseif(kx > (numl - 2))
%         kx = numl - 1;
%     end
%     if(ky < 0)
%         ky = -1;
%     elseif(ky > (numl - 2))
%         ky = numl - 1;
%     end
%     rxmin = tl * (2 * kx + 1) + tsdl;
%     rxmax = tl * (2 * kx + 3) + tsdl;
%     rymin = tl * (2 * ky + 1) + tsdl;
%     rymax = tl * (2 * ky + 3) + tsdl;
%     a=[];
%     a(1) = circleSDF(x,[rxmin,rymin],Rl);
%     a(2) = circleSDF(x,[rxmax,rymin],Rl);
%     a(3) = circleSDF(x,[rxmin,rymax],Rl);
%     a(4) = circleSDF(x,[rxmax,rymax],Rl);
%     a(5) = t;
%     a(6) = t;
%     a(7) = t;
%     a(8) = t;
%     %   边缘的情况
%     if(kx < 0)% rxmin的圆不能要了（不存在）
%         a(1) = t;
%         a(3) = t;
%         rxmin = -tr - tsdr;
%         kky = floor(((x(2) - tsdr)/tr-1)/2);
%         rymin = tr * (2 * kky + 1) + tsdr;
%         rymax = tr * (2 * kky + 3) + tsdr;
%         a(5) = circleSDF(x,[rxmin,rymin],Rr);
%         a(6) = circleSDF(x,[rxmin,rymax],Rr);
%         if(kky < 0)
%             a(5) = t;
%         elseif(kky > (numr - 2))
%             a(6) = t;
%         end
%     elseif(kx > (numl - 2))% rxmax的圆不能要了（不存在）
%         a(2) = t;
%         a(4) = t;
%         rxmax = t + tr + tsdr;
%         kky = floor(((x(2) - tsdr)/tr-1)/2);
%         rymin = tr * (2 * kky + 1) + tsdr;
%         rymax = tr * (2 * kky + 3) + tsdr;
%         a(5) = circleSDF(x,[rxmax,rymin],Rr);
%         a(6) = circleSDF(x,[rxmax,rymax],Rr);
%         if(kky < 0)
%             a(5) = t;
%         elseif(kky > (numr - 2))
%             a(6) = t;
%         end
%     end
%     if(ky < 0)% rymin的圆不能要了（不存在）
%         a(1) = t;
%         a(2) = t;
%         rymin = -tr - tsdr;
%         kkx = floor(((x(1) - tsdr)/tr-1)/2);
%         rxmin = tr * (2 * kkx + 1) + tsdr;
%         rxmax = tr * (2 * kkx + 3) + tsdr;
%         a(7) = circleSDF(x,[rxmin,rymin],Rr);
%         a(8) = circleSDF(x,[rxmax,rymin],Rr);
%         if(kkx < 0)
%             a(7) = t;
%         elseif(kkx > (numr - 2))
%             a(8) = t;
%         end
%     elseif(ky > (numl - 2))% rymax的圆不能要了（不存在）
%         a(3) = t;
%         a(4) = t;
%         rymax = t + tr + tsdr;
%         kkx = floor(((x(1) - tsdr)/tr-1)/2);
%         rxmin = tr * (2 * kkx + 1) + tsdr;
%         rxmax = tr * (2 * kkx + 3) + tsdr;
%         a(7) = circleSDF(x,[rxmin,rymax],Rr);
%         a(8) = circleSDF(x,[rxmax,rymax],Rr);
%         if(kkx < 0)
%             a(7) = t;
%         elseif(kkx > (numr - 2))
%             a(8) = t;
%         end
%     end
%     sd = min(a);
% end
%%  bulk-porous media体系,bulk为一个圆，粒子以某个频率向中间发射，中间是多孔介质的体系，也算作方形阵列
%  R = [r,Rp,nump],r为bulk圆半径,圆心在原点,Rp为阵列的孔喉比,nump为圆阵列每行个数，圆的半径为1
% if(length(R) ~= 3)
%     error('孔隙参数设置错误');
% end
% r = R(1);
% Rp = R(2);
% nump = R(3);
% sq2 = 2^0.5;
% t = nump * Rp;
% if(t * sq2 > r)
%     error('bulk圆比圆阵列还小');
% end
% tsd = (x(1)^2 + x(2)^2)^0.5;
% if(tsd > r + 0.1)
%     sd = nan;
% elseif(tsd > t * sq2)
%     sd = tsd - t * sq2 + Rp - 1;
% else
%     x = x + t;
%     kx = floor((x(1)/Rp-1)/2);
%     ky = floor((x(2)/Rp-1)/2);
%     if(kx < 0)
%         kx = -1;
%     elseif(kx > (nump - 2))
%         kx = nump - 1;
%     end
%     if(ky < 0)
%         ky = -1;
%     elseif(ky > (nump - 2))
%         ky = nump - 1;
%     end
%     rxmin = Rp * (2 * kx + 1);
%     rxmax = Rp * (2 * kx + 3);
%     rymin = Rp * (2 * ky + 1);
%     rymax = Rp * (2 * ky + 3);
% 
%     a=[];
%     a(1) = circleSDF(x,[rxmin,rymin],1);
%     a(2) = circleSDF(x,[rxmax,rymin],1);
%     a(3) = circleSDF(x,[rxmin,rymax],1);
%     a(4) = circleSDF(x,[rxmax,rymax],1);
%     if(kx < 0)% rxmin的圆不能要了（不存在）
%         a(1) = t;
%         a(3) = t;
%     elseif(kx > (nump - 2))% rxmax的圆不能要了（不存在）
%         a(2) = t;
%         a(4) = t;
%     end
%     if(ky < 0)% rymin的圆不能要了（不存在）
%         a(1) = t;
%         a(2) = t;
%     elseif(ky > (nump - 2))% rymax的圆不能要了（不存在）
%         a(3) = t;
%         a(4) = t;
%     end
%     sd = min(a);
% end

%%  圆形空腔
% if(length(R) ~= 1)
%     error('孔隙参数设置错误');
% end
% r = R(1);
% tsd = (x(1)^2 + x(2)^2)^0.5;
% if(tsd > r + 0.1)
%     sd = nan;
% else
%     sd = tsd + 0.1;
% end
%%  圆形+孔隙,R=[大圆半径，小圆半径，孔隙角度],圆弧厚度1e-3
% if(length(R) ~= 3)
%     error('孔隙参数设置错误');
% end
% if(R(3)/2 > pi/6)
%     error('角度设置不对');
% end
% r = R(1);
% rr = R(2);
% TH = R(3)/2;
% flag = 0;
% tsd = (x(1)^2 + x(2)^2)^0.5;
% th = atan2(x(2), x(1));
% if((abs(th - pi/2) < TH) || (abs(th + pi/6) < TH) || (abs(th + 5*pi/6) < TH))
%     flag = 1;
% end
% if(abs(th + pi/6) < TH)
%     x = x * [-1/2 3^0.5/2;-3^0.5/2 -1/2];
% end
% if(abs(th + 5*pi/6) < TH)
%     x = x * [-1/2 -3^0.5/2;3^0.5/2 -1/2];
% end
% k = tan(pi/2-TH);
% s = min(abs(k*x(1)+x(2))/(1+k^2),abs(k*x(1)-x(2))/(1+k^2));
% x0 = (x(1)/k+x(2))/(k+1/k);
% y0 = k * x0;
% k = -k;
% x01 = (x(1)/k+x(2))/(k+1/k);
% y01 = k * x0;
% ttsd = (x0 ^ 2 + y0 ^ 2)^0.5;
% ttsd1 = (x01 ^ 2+y01 ^ 2)^0.5;
% ttsd = max(ttsd,ttsd1);
% if(tsd > r + 0.1)
%     sd = nan;
% elseif(flag == 0)
%     if(tsd > rr)
%         sd = tsd - rr;
%     elseif(tsd < rr - 1e-3)
%         sd = rr - 1e-3 - tsd;
%     else
%         sd = min(rr - tsd, tsd - rr + 1e-3);
%     end
% else
%     if(ttsd > rr)
%         sd = ttsd - rr;
%         sd = (sd^2 + s^2)^0.5;
%     elseif(ttsd < rr - 1e-3)
%         sd = rr - 1e-3 - ttsd;
%         sd = (sd^2 + s^2)^0.5;
%     else
%         sd = s;
%     end
% end
%%  圆阵列随机场
%   cir_1:420*420多孔介质，1300个圆，500的半径
%   cir_2:四个圆组成一个孔，R=1.1
% cir = load('data\cir.mat');
% cir = cir.cir_2;
% if(length(R) ~= 1)
%     error('孔隙参数设置错误');
% end
% r = R(1);
% len = size(cir);len = len(1);
% s = zeros(1,len);
% tsd = (x(1)^2 + x(2)^2)^0.5;
% if(tsd > r + 0.1)
%     sd = nan;
% elseif(tsd > 220*2^0.5)
%     sd = tsd - 220*2^0.5 + 0.1;
% else
%     for i = 1:len
%         s(i) = circleSDF(x, [cir(i,1), cir(i,2)], cir(i,3));
%     end
%     sd = min(s);
% end
%%  三角界面
% s = [84 * 3^0.5, 1];
% if(length(R) ~= 1)
%     error('孔隙参数设置错误');
% end
% r = R(1);
% tsd = (x(1)^2 + x(2)^2)^0.5;
% if(tsd > r + 0.1)
%     sd = nan;
% else
%     sd = min([boxSDF(x, [0, -105], s, 0), boxSDF(x, [-52.5 * 3^0.5, 52.5], s, pi/3), boxSDF(x, [52.5 * 3^0.5, 52.5], s, -pi/3)]);
% end    
%%  交叉圆
% if(length(R) ~= 1)
%     error('孔隙参数设置错误');
% end
% r = R(1);
% tsd = (x(1)^2 + x(2)^2)^0.5;
% if(tsd > r + 0.1)
%     sd = nan;
% else
%     sd = min(arcSDF(x, [-150, 0], -50/180*pi, 50/180*pi, 200, 205), arcSDF(x, [150, 0], pi-50/180*pi, pi+50/180*pi, 200, 205));
% end
%%  四个小正方形
% if(length(R) ~= 1)
%     error('孔隙参数设置错误');
% end
% r = R(1);
% tsd = (x(1)^2 + x(2)^2)^0.5;
% if(tsd > r + 0.1)
%     sd = nan;
% else
%     sd = min([boxSDF(x, [150, 150], [30, 30], 0), boxSDF(x, [-150, 150], [30, 30], 0), boxSDF(x, [150, -150],...
%          [30, 30], 0), boxSDF(x, [-150, -150], [30, 30], 0)]);
% end
%%  正方形阵列
%  R = [r,Rp,nump],r为bulk圆半径,圆心在原点,Rp为阵列的(正方形边长+间距)/边长比,nump为正方形阵列每行个数，正方形的半边长为1
% if(length(R) ~= 3)
%     error('孔隙参数设置错误');
% end
% s = [1 1];
% r = R(1);
% Rp = R(2);
% nump = R(3);
% sq2 = 2^0.5;
% t = nump * Rp;
% if(t * sq2 > r)
%     error('bulk圆比圆阵列还小');
% end
% tsd = (x(1)^2 + x(2)^2)^0.5;
% if(tsd > r + 0.1)
%     sd = nan;
% elseif(tsd > t * sq2)
%     sd = tsd - t * sq2 + Rp - 1;
% else
%     x = x + t;
%     kx = floor((x(1)/Rp-1)/2);
%     ky = floor((x(2)/Rp-1)/2);
%     if(kx < 0)
%         kx = -1;
%     elseif(kx > (nump - 2))
%         kx = nump - 1;
%     end
%     if(ky < 0)
%         ky = -1;
%     elseif(ky > (nump - 2))
%         ky = nump - 1;
%     end
%     rxmin = Rp * (2 * kx + 1);
%     rxmax = Rp * (2 * kx + 3);
%     rymin = Rp * (2 * ky + 1);
%     rymax = Rp * (2 * ky + 3);
% 
%     a=[];
%     a(1) = boxSDF(x,[rxmin,rymin],s,0);
%     a(2) = boxSDF(x,[rxmax,rymin],s,0);
%     a(3) = boxSDF(x,[rxmin,rymax],s,0);
%     a(4) = boxSDF(x,[rxmax,rymax],s,0);
%     if(kx < 0)% rxmin的圆不能要了（不存在）
%         a(1) = t;
%         a(3) = t;
%     elseif(kx > (nump - 2))% rxmax的圆不能要了（不存在）
%         a(2) = t;
%         a(4) = t;
%     end
%     if(ky < 0)% rymin的圆不能要了（不存在）
%         a(1) = t;
%         a(2) = t;
%     elseif(ky > (nump - 2))% rymax的圆不能要了（不存在）
%         a(3) = t;
%         a(4) = t;
%     end
%     sd = min(a);
% end
%%  三凹弧
% if(length(R) ~= 1)
%     error('孔隙参数设置错误');
% end
% r = R(1);
% tsd = (x(1)^2 + x(2)^2)^0.5;
% if(tsd > r + 0.1)
%     sd = nan;
% else
%     sd = min([arcSDF(x, [0 -150], pi/2-50/180*pi, pi/2+50/180*pi, 120, 125), ...
%         arcSDF(x, [150/2*3^0.5, 150/2], 2*pi/3+pi/2-50/180*pi, 2*pi/3+pi/2+50/180*pi, 120, 125) ...
%         arcSDF(x, [-150/2*3^0.5, 150/2], -2*pi/3+pi/2-50/180*pi, -2*pi/3+pi/2+50/180*pi, 120, 125)]);
% end
%%  三凸弧
% if(length(R) ~= 1)
%     error('孔隙参数设置错误');
% end
% r = R(1);
% tsd = (x(1)^2 + x(2)^2)^0.5;
% if(tsd > r + 0.1)
%     sd = nan;
% else
%     sd = min([arcSDF(x, [0 -200], pi/2+30/180*pi, pi/2-30/180*pi, 150, 155), ...
%         arcSDF(x, [200/2*3^0.5, 200/2], 2*pi/3+pi/2+30/180*pi, 2*pi/3+pi/2-30/180*pi, 150, 155) ...
%         arcSDF(x, [-200/2*3^0.5, 200/2], -2*pi/3+pi/2+30/180*pi, -2*pi/3+pi/2-30/180*pi, 150, 155)]);
% end
%%  随机矩形多孔介质
% cir = load('cube.mat');
% cir = cir.cube_1;
% if(length(R) ~= 1)
%     error('孔隙参数设置错误');
% end
% r = R(1);
% len = size(cir);len = len(1);
% s = zeros(1,len);
% tsd = (x(1)^2 + x(2)^2)^0.5;
% if(tsd > r + 0.1)
%     sd = nan;
% elseif(tsd > 220*2^0.5)
%     sd = tsd - 220*2^0.5 + 0.1;
% else
%     for i = 1:len
%         s(i) = boxSDF(x, [cir(i,1), cir(i,2)], [cir(i,3), cir(i,4)], cir(i,5));
%     end
%     sd = min(s);
% end
%%  引用scene_pix
ppm = p.Results.pm;ppm = ppm{1};
if(isequal(ppm,0))
    error('pm没给');
end
if(length(R) ~= 1)
    error('孔隙参数设置错误');
end
r = R(1);
t = 205;
sq2 = 2^0.5;
if(t * sq2 > r)
    error('bulk圆比圆阵列还小');
end
tsd = (x(1)^2 + x(2)^2)^0.5;
tssd = boxSDF(x, [0 0], [t t], 0);
if(tsd > r + 0.1)
    sd = nan;
elseif(tssd > 0)
    sd = tssd + 1;
else
    sd = scene_psi(x(1), x(2), ppm.sdf, [-t t], [-t t], ppm.maxx, ppm.maxy);
end














end
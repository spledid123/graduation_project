%   判断"入口效应"的程序
%   一根直管，等距离分成孔，统计单孔碰撞次数,统计分布与孔长度的关系（孔越长，离指数越远）
%   统计跳入某个孔后出来回到原方向的概率判断"入口效应"的大小及与在孔内碰撞次数的关系；与孔的长度的关系
%   概率逐渐降为0.5，"入口效应"消失

%   输入
ts = 0.5;                   %  孔的长度
num = 1e6;                  %  碰撞次数

wid = 1;                    %  直管宽度
Dtp = [];                   %  单孔碰撞次数
fun_cosin=@(x)asin(2.*x-1); %  余弦定律
fun_pore=@(x)floor(x/ts);   %  判断孔编号
A=zeros(1,num+1);           %  碰撞位置
Dk = [];                    %  每次碰撞的孔编号
A(1)=0.9;                   %  初始化，定起点
Dk(1) =fun_pore(A(1));
Dtp(1)=1;
%   模拟
j=1;
for i=1:num
    A(i+1)=A(i)+wid*tan(fun_cosin(rand()));
    if(fun_pore(A(i+1))~=k)%    越孔
        j = j + 1;
        Dtp(j) = 1;
        Dk(j) = fun_pore(A(i+1)/ts);
    else                    %   不越孔
        Dtp(j) = Dtp(j) + 1;
    end
end
%   统计，画频率图
Dr1 = zeros(1,max(Dtp));%   统计跳入孔内在孔内碰撞k次后回到原方向的频数
Dr2 = zeros(1,max(Dtp));%   统计跳入孔内在孔内碰撞k次后没有回到起点方向的频数
for i=2:length(Dk)-1
    if((Dk(i-1)<Dk(i))==(Dk(i+1)<Dk(i)))%   从i-1表示的孔->i表示的孔->i+1表示的孔；i-1表示的孔与i+1表示的孔在同一边
        Dr1(Dtp(i))=Dr1(Dtp(i))+1;
    else
        Dr2(Dtp(i))=Dr2(Dtp(i))+1;
    end
end
Dr = Dr1 + Dr2;
Pr1 = Dr1./Dr;
figure;
hold on;
box on;
plot(Pr1);% 跳入孔内在孔内碰撞k次后回到原方向的概率
title('出发前往是同一边的概率','fontsize',14);
xlabel('单孔碰撞次数','FontSize',14);
ylabel('回到同一边的概率','FontSize',14);
%%   单孔碰撞次数分布
[xx,yy] = plot_distribution(Dtp,1,10);
figure;
hold on;
box on;
plot(xx,yy);
set(gca,'yscale','log');
set(gca,'ylim',[max(Ptt)/1e2 max(Ptt)]);
title('单孔停留次数概率','FontSize',14);
xlabel('单孔停留次数','FontSize',14);
ylabel('概率','FontSize',14);
%   画出某次跨越边界的轨迹
%   需要补充，修改
A = readtable("data\partical_longpath_cir\rxT_circle_onepartical_cubecirR_1.1.txt");
rx = table2array(A(1:end,2));
ry = table2array(A(1:end,3));
j = 1;
kAA = [];
hold on;
% 筛选
for i=1:length(kTT)
    kAA(j)=i;
    j=j+1;
end
ki=[];
for i=1:length(kAA)
    ki(i)=kA(kAA(i));
    i=i+1;
end
f=19;
kxx = floor((rx(ki(f))/R(1)-1)/2);
kyy = floor((ry(ki(f))/R(2)-1)/2);
kx=[kxx-2 kxx-1 kxx kxx+1 kxx+2];
ky=[kyy-2 kyy-1 kyy kyy+1 kyy+2];
numx = 5;
numy = 5;
for i=1:numx
    x0(i)=R(1)*(2*kx(i)+1);
end
for i=1:numy
    y0(i)=R(2)*(2*ky(i)+1);
end
hold on
for i=1:numx
    for j=1:numy
         plot_cir(x0(i),y0(j),1);
    end
end
plot(x0(3),y0(3),'or');
axis equal
rxx = rx(ki(f)-50:ki(f)+2);
ryy = ry(ki(f)-50:ki(f)+2);
plot(rxx,ryy);



%   生成随机数
%   随机数满足3D/2D玻尔兹曼分布
%   储存在相应文件中
%   需要更改Boltzmann.m文件更改分布
filename = 'data\v_10000\Blotzmann_water_T_40.txt';
sumv = 1e4;
v = 0;
varNames = {'v'};
for i = 1:sumv
    v = Boltzmann(40);
    writetable(table(v,'VariableNames',varNames),filename,'WriteMode','Append', 'WriteVariableNames',false,'WriteRowNames',true);
end

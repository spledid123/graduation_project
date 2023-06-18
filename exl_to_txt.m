%%   之前保存文件一直使用xlsx格式，但这样读取写入效率很慢。
%   批量读取文件夹下的所有xlsx文件，将其改为txt文件
%   使用readtable和writetable
a = dir('**/*.xlsx');% 读取文件夹里所有xlsx文件(包括子文件夹)
len = length(a);
parfor i = 1:len
    filename = strcat(a(i).folder,'\',a(i).name);
    A = readtable(filename,'ReadVariableNames',true,'NumHeaderLines',1);%  之前的变量名字一直在第2行
    old = '.xlsx';
    new = '.txt';
    newfilename = strrep(filename,old,new);%    替换字符串
    writetable(A,newfilename,'WriteMode','overwrite','WriteVariableNames',true);    
    if(mod(i,100) == 0)
        fprintf('文件转化进度为%f%%\n',i/len*100);
    end
    delete(filename);%  删除文件
end
%%  删除文件夹里不需要的文件
a = dir('data/bulk_pore_pm13_R_500_N_1000000_dT_1_2000/*.txt');% 读取文件夹里所有txt文件(包括子文件夹)
len = length(a);
for i = 1:len
    t = a(i).name;
    pat = digitsPattern;
    t = extract(t,pat);
    t = t{1};
    t = str2double(t);
    if(mod(t,10)~=0)
        filename = strcat(a(i).folder,'\',a(i).name);
        delete(filename);
    end
end



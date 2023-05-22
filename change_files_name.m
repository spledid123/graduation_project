%   文件名批处理
filenamesta="data\path_2000_six_R_1.1_a_0.25(2)\rxT_circle_onepartical_sixcir_";
filenameend = '_.txt';
for m=1:500
    filenamemid = int2str(m);
    filename = strcat(filenamesta,filenamemid,filenameend);
    filenamemid = int2str(m+500);
    newname = strcat(filenamemid,filenameend);
    copyfile([filename],['data\path_2000_six_R_1.1_a_0.25\rxT_circle_onepartical_sixcir_' newname]);% 可使用movefile作剪切
end
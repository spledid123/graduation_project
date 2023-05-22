function [] = plot_spe(t, x)
%   画出频谱图，找周期
%   函数为f(t)=x
%   t为等距采样，为等差数列
len = length(t);
su = (t(1)+t(len))*len/2;
if(abs(su-sum(t))>1e-3)
    error('xx不是等差数列');
end
if(length(x) ~= len)
    error('xx,yy长度不同');
end
newx = x - mean(x);
y = fft(newx);
y = abs(y);
figure;
plot((1:floor(len/2))/(t(end) - t(1)), y(1:floor(len/2)));
xlabel('Frequency');
ylabel('Magnitude');

end
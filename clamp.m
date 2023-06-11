function sd = clamp(value, minvalue, maxvalue)
% 如果参数位于最小数值和最大数值之间的数值范围内，则该函数将返回参数值。
% 如果参数大于范围，该函数将返回最大数值。
% 如果参数小于范围，该函数将返回最小数值。
if(minvalue > maxvalue)
    error('最大值小于最小值');
end
if(value < minvalue)
    sd = minvalue;
elseif(value > maxvalue)
    sd = maxvalue;
else
    sd = value;
end
end
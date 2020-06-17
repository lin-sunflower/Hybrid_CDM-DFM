function [b1,b2] = decodeGRAY(symbol)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
for c=1:1:length(symbol)
    if symbol(c)==1
        b1(c)=0;b2(c)=0;
    end
        if symbol(c)==2
        b1(c)=0;b2(c)=1;
        end
        if symbol(c)==3
        b1(c)=1;b2(c)=1;
        end
        if symbol(c)==4
        b1(c)=1;b2(c)=0;
    end
end
end


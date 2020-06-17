function [symbol] = Graymap(seq,entropy)
i=1;
for c=1:2:floor(length(seq)/2)
    if seq(c)==-1 && seq(c+1)==-1
        symbol(i)=0;
        i=i+1;
    end
    if seq(c)==-1 && seq(c+1)==1
        symbol(i)=1;
        i=i+1;
    end
    if seq(c)==1 && seq(c+1)==-1
        symbol(i)=2;
        i=i+1;
    end
    if seq(c)==1 && seq(c+1)==1
        symbol(i)=3;
        i=i+1;
    end    
end
symbol=symbol-1.5; %dc block.
end

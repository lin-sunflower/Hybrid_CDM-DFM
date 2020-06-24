function BER = theore_ber(inputsymbol,outputsymbol,start)
Length=min([length(inputsymbol) length(outputsymbol)]);
inputsymbol=inputsymbol(start:Length);
outputsymbol=outputsymbol(start:Length);
power=sum(outputsymbol.^2)/(length(outputsymbol));
outputsymbol=outputsymbol/sqrt(power);
for i=1:1:4
STD(i)=std(outputsymbol(find(inputsymbol==(i-1-1.5))));
Level(i)=mean(outputsymbol(find(inputsymbol==(i-1-1.5))));
end
for k=2:4
    I_th(k)=[Level(k)+Level(k-1)]/2;
end
I_th(1)=-inf;
I_th(5)=inf;

for x=1:4       %transmit x
    for y=1:4   %received y
        %probability of receiving y while x was transmitted
        Prob(x,y)=0.5*erfc([I_th(y)-Level(x)]/[sqrt(2)*STD(x)])-0.5*erfc([I_th(y+1)-Level(x)]/[sqrt(2)*STD(x)]);
        Prob(x,y)=1/4*Prob(x,y);
    end
end
for k=1:4
   Prob(k,k)=0;   
end
BER=1/log2(4)*sum(Prob(:));

end


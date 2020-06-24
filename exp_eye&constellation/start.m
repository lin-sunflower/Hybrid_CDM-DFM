function [start] = start(symbolrate, samplerate, tx_data, rxdata)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
label=tx_data;
data=rxdata;
sample=symbolrate;  
LL=200000;
k=30000;
T=length(label);

tem=(max(label)-min(label))/(max(data)-min(data));
data=data*tem;

M1=mean(label);
M2=mean(data);
data=data+M1-M2;

Temp=ceil(length(data)/(T*160/sample));

for i=1:Temp
    label((i-1)*T+1:i*T)=label(1:T);
end  %构造足够的类别数量

n_l=1:length(label);  
m_l=samplerate/sample;   
s_l=1/m_l;         
xi_l=1:s_l:length(label);
label1= interp1(n_l,label,xi_l,'linear');

for i=1:LL
    sdata=data(i:k+i-1);
    Error=sdata-label1(1:length(sdata));
    error(i)=var(Error);
end
[vlu1 ind1]=sort(error);

start=ind1(1)
figure;plot(data(start:32768+start-1));hold on ;plot(label1(1:32768));
end


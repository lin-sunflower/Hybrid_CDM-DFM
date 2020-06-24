clear;clc;close all;
%------- this file for data generation and processing for experimental use
M=4;%%%%%%%%%% CDM lane
N=2;%%%%%%%%%% DFM lane
upsampling=3;
%%%%%------- DFM shaping ----------
a_cap=0.1;                  %% Upsampling
cap_length=50;
Fd_cap=1/3;
Fs_cap=1;
w_cap=rcosine(Fd_cap, Fs_cap, 'fir/sqrt', a_cap, cap_length);
x_cap=-cap_length:Fd_cap/Fs_cap:cap_length;
f1=w_cap.*sin(1.1*pi*x_cap);
f2=w_cap.*cos(1.1*pi*x_cap);
f=[f1;f2];
k_cap=zeros(1,Fs_cap/Fd_cap-1);k_cap=[1 k_cap]';
%------ information symbol --------
L=13;
binary=idinput(2^L -1, 'prbs'); 
symbol=Graymap(binary,2);
%-------   CDM-DFM modulation  ------
Walshcode=hadamard(M);
for c1=1:1:N
for c2=1:1:M
    code_symbol_I(c2,:)=ds_mod(Walshcode(c2,:)',symbol((c1-1)*M+c2:M*N:end));
end
Total_code_I=sum(code_symbol_I);
up_Symbol=kron(Total_code_I,k_cap');
filter_data(c1,:)=(-1)^c1*filter(f(c1,:),1,up_Symbol);
end
Tx_data=sum(filter_data);
% Tx_data=resample(Tx_data,92,20);
% %------  Power normalization  -----
Tx_data=Tx_data-mean(Tx_data);
power=sum(Tx_data.^2)/(length(Tx_data));
Tx_data=Tx_data/sqrt(power);
Tx_data=Tx_data-min(Tx_data);
Tx_data=sqrt(Tx_data);  
%-----   Rx processing on one frame ------
Rx_symbol=Tx_data.^2;
% save cdm4dfm2.mat Tx_data
%-----   CDM decoding  ------
for c1=1:1:N
defilter_data=filter(f(c1,:),1,Rx_symbol);
defilter_S=defilter_data(length(w_cap):Fs_cap/Fd_cap:end);
defilter_symbol(c1,:)=defilter_S(1:floor(length(defilter_S)/4)*4);
for c2=1:1:M
    decode_symbol((c1-1)*M+c2,:)= ds_demod(Walshcode(c2,:)',defilter_symbol(c1,:)); 
end
end
%-----   scatter plots  ------
lane=1;
scatter(decode_symbol(lane,:),decode_symbol(lane+M,:))
for c=1:1:M*N
    ber(c)=theore_ber(symbol(c:M*N:end),decode_symbol(c:M*N:end),1);
end
BER_total=sum(ber)/(M*N);

BER_total=BER_total

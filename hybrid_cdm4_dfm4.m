clear;clc;
% figure;
M=4;%%%%%%%%%% CDM lane
N=4;%%%%%%%%%% DFM lane
upsampling=5;
%------ channel rate setup --------
lw=20e6;     %% linewidth
Fd=56e9;     %%% symbol rate in the physical layer
Fs=Fd*30/upsampling;     %% sampling rate
BW=20e9;
H=gaussfir(1/(Fs)*BW,50,30);
%%%%%------- DFM shaping ----------
a_cap=0.1;                  %% Upsampling
cap_length=20;
Fd_cap=1/upsampling;
Fs_cap=1;
w_cap=rcosine(Fd_cap, Fs_cap, 'fir/sqrt', a_cap, cap_length);
x_cap=-cap_length:Fd_cap/Fs_cap:cap_length;
f1_i=w_cap.*sin((2*1-1+0.1)*pi*x_cap);
f1_q=w_cap.*cos((2*1-1+0.1)*pi*x_cap);
f2_i=w_cap.*sin((2*2-1+0.25)*pi*x_cap);
f2_q=w_cap.*cos((2*2-1+0.25)*pi*x_cap);
f=[f1_i;f1_q;f2_i;f2_q];
k_cap=zeros(1,Fs_cap/Fd_cap-1);k_cap=[1 k_cap]';
%-------      DA       ------
a=0.5;                     
DA_length=20;
w_DA=rcosine(Fd, Fs, 'fir/sqrt', a, DA_length);
k=zeros(1,Fs/Fd-1);k=[1 k];
%------ information symbol --------
L=21;
binary=idinput(2^L -1, 'prbs'); 
symbol=Graymap(binary,2);
%------ MPI configuration setup --------
% Length=25;
Length=10;
for alpha=0:0.1:0.2
for j=1:1:20
SNR(j)=j;   
snr=10^(SNR(j)/10);
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
Tx_symbol=sum(filter_data);
%-------      DA       ------
up_Symbol=kron(Tx_symbol,k);
Tx_data=filter(w_DA,1,up_Symbol);
%------  Transmission  -----
Tx_data=Tx_data-mean(Tx_data);
power=sum(Tx_data.^2)/(length(Tx_data));
Tx_data=Tx_data/sqrt(power);
Tx_data=Tx_data-min(Tx_data);
Tx_data=sqrt(Tx_data);     
%   MPI and gaussian noise
c_light=3e8;
n_1=1.45;
Td=round(2*Length*n_1/c_light*Fs);
wc=2*pi*192e12;
N_phase=0 + sqrt(2*pi*lw*Td/Fs).*randn(length(Tx_data),1);
N_inten=0 + sqrt(1/snr/2).*randn(length(Tx_data),1);
for t=1:1:Td
    Rx_data(t)=Tx_data(t)^2+N_inten(t);
end
for t=Td+1:1:length(Tx_data)
   Rx_data(t)=Tx_data(t)^2+alpha^2*Tx_data(t-Td)^2+2*alpha*Tx_data(t)*Tx_data(t-Td)*cos(wc*Td+N_phase(t))+N_inten(t);
end
%--------------  bandwidth limitation  ----------------
K=zeros(1,29);K=[1 K];
Rx_data=kron(Rx_data,K);
Rx_data=conv(Rx_data,H);
Rx_data=Rx_data((length(H)-1)/2+1:30:end-length(H)/2);  %%% 20 is the sample resolution of channel
%-------      AD       ------
Rx_data=filter(w_DA,1,Rx_data);
Rx_symbol=Rx_data(length(w_DA):Fs/Fd:end);
%-----   CDM decoding  ------
for c1=1:1:N
defilter_data=filter(f(c1,:),1,Rx_symbol);
defilter_S=defilter_data(length(w_cap):Fs_cap/Fd_cap:end);
defilter_symbol(c1,:)=defilter_S(1:floor(length(defilter_S)/4)*4);
for c2=1:1:M
    decode_symbol((c1-1)*M+c2,:)= ds_demod(Walshcode(c2,:)',defilter_symbol(c1,:)); 
end
end
Start=ceil(Td/(Fs/Fd));
for c=1:1:M*N
    ber(c)=theore_ber(symbol(c:M*N:end),decode_symbol(c:M*N:end),Start);
end
BER_total(j)=sum(ber)/(M*N);
end
% figure;
hold on; plot(SNR,log10(BER_total),'-')
xlabel('SNR(dB)');ylabel('log_{10}(BER)');
grid on;
end

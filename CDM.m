clear;clc;
% close all;
% figure;
hold on;
M=8;%%%%%%%%%% CDM lane 8 and 16
%------ channel rate setup --------
lw=20e6;     %% linewidth
Fd=56e9;     %%% symbol rate in the physical layer
Fs=Fd*30;     %% sampling rate
BW=20e9;
H=gaussfir(1/(Fs)*BW,50,30);
%------ DFM or CAP setup --------
entropy=4;   %% qam-16
upsampling_second=2;
band=1;
information_rate=Fd/upsampling_second*band*entropy/1e9;
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
%-------   CDM coding  ------
Walshcode=hadamard(M);
for c=1:1:M
    code_symbol(c,:)=ds_mod(Walshcode(c,:)',symbol(c:M:end));
end
Tx_symbol=sum(code_symbol);
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
Td=round(Length*n_1/c_light*Fs);
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
Rx_symbol=Rx_symbol(1:floor(length(Rx_symbol)/M)*M);
for c=1:1:M
    decode_symbol(c,:)=ds_demod(Walshcode(c,:)',Rx_symbol);
end
Start=ceil(Td/(Fs/Fd));
for c=1:1:M
    ber(c)=theore_ber(symbol(c:M:end),decode_symbol(c,:),Start)
end
BER_total(j)=sum(ber)/M;
end
hold on;plot(SNR,log10(BER),'-*')
xlabel('SNR(dB)');ylabel('log_{10}BER');
grid on;
end

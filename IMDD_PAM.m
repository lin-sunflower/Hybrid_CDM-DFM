clear;clc;
% close all;
%------ channel sampling rate setup --------
lw=20e6;     %% linewidth
Fd=56e9;     %%% symbol rate in the physical layer
Fs=Fd*30;     %% sampling rate
BW=20e9;
H=gaussfir(1/(Fs)*BW,50,30);
% figure;plot(H)
%------ DFM or CAP setup --------
entropy=4;   %% qam-16
upsampling_second=2;
band=1;
information_rate=Fd/upsampling_second*band*entropy/1e9;
%%%%%------- DA ----------
a=0.5;                
Fir_length=20;
w=rcosine(Fd, Fs, 'fir/sqrt', a, Fir_length);
k=zeros(1,Fs/Fd-1);k=[1 k]';
%------ MPI configuration setup --------
% Length=25;
Length=10; %%%%% 30 graduality
% subplot(1,2,1);
for alpha=0:0.1:0
%--------------------------------
for j=1:1:20
SNR(j)=j;   %% Intensity noise in dB
snr=10^(SNR(j)/10);
L=21;
Binary=idinput(2^L -1, 'prbs'); 
Symbol=Graymap(Binary,2);
%%%%%------- Shaping with SRRC waveform----------
up_Symbol=kron(Symbol',k);
Tx_data=filter(w,1,up_Symbol);
%%%%%------- normalized power ----------
Tx_data=Tx_data-mean(Tx_data);
power=sum(Tx_data.^2)/(length(Tx_data));
Tx_data=Tx_data/sqrt(power);
Tx_data=Tx_data-min(Tx_data);
Tx_data=sqrt(Tx_data);     %% For square-law detection (MZ)
%%%%%------- MPI by 2-times reflection----------
c=3e8;
n_1=1.45;
Td=round(2*Length*n_1/c*Fs);
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
%%%%%------- sampling and decoding----------
Rx_data=filter(w,1,Rx_data);
Rx_symbol=Rx_data(length(w):Fs/Fd:end);
% eyediagram(Rx_data(length(w)-1:50000),Fs/Fd,1,1);       %%% eyediagram has to shit to 'time 0' by mimus 1.
% eyediagram(Rx_symbol(1:10000),2,1,1);
%%%%%------- BER----------
Start=ceil(Td/(Fs/Fd));
BER(j)=theore_ber(Symbol,Rx_symbol,Start);
end
hold on;plot(SNR,log10(BER),'-*')
xlabel('SNR(dB)');ylabel('log_{10}BER');
grid on;
end
clear;clc;close all;
symbolrate=84;
samplerate=160;
M=4;%%%%%%%%%% CDM lane
% M=8;%%%%%%%%%% CDM lane
N=2;%%%%%%%%%% DFM lane
%---- this file for data receiving and processing for experimental use
load('cdm4dfm2.mat');
% load('cdm8dfm2.mat');
length_frame = length(Tx_data);
%---- this paragraph for sample datas on DSO
% addpath('scope')
% wav = superman(160e9, 10e6)';
%------------------------
load rx_dso4(a)0.mat
% save rx_dso4_BI(a).mat wav
% save rx_dso4(b).mat wav
% save rx_dso4_BI(b).mat wav

% save rx_dso8(a).mat wav
% save rx_dso8_BI(a).mat wav
% save rx_dso8(b).mat wav
% save rx_dso8_BI(b).mat wav
%------------------------
L=13;
binary=idinput(2^L -1, 'prbs');
symbol=Graymap(binary,2);
%------------------------
point=start(symbolrate, samplerate, Tx_data.^2, wav);

wav=wav(point:length(wav)/5);
wav = downsample_data(wav, symbolrate, samplerate, 0.1);
% figure;plot(wav);hold on;plot(Tx_data)

wav = LMS_FFE(Tx_data.^2, wav, 121, 0.0001);
% figure;plot(wav);hold on;plot(Tx_data.^2)
%-------------- refer to Tx settings  ---------------------
Walshcode=hadamard(M);
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
% figure;plot(f1);hold on;plot(f2)
% figure;plot(10*log10(abs(fft(f1))));
% hold on;plot(10*log10(abs(fft(f2))));
%-----   Rx processing for each frame ------

num_frame = floor(length(wav)/length(Tx_data));
for c_frame=1:1:num_frame
    Rx_symbol=wav(1+(c_frame-1)*length_frame:c_frame*length_frame);
    
    %-----   CDM decoding  ------
    for c1=1:1:N
        defilter_data=filter(f(c1,:),1,Rx_symbol);
        defilter_S=defilter_data(length(w_cap):Fs_cap/Fd_cap:end);
        defilter_symbol(c1,:)=defilter_S(1:floor(length(defilter_S)/M)*M);
        for c2=1:1:M
            decode_symbol((c1-1)*M+c2,c_frame,:)= ds_demod(Walshcode(c2,:)',defilter_symbol(c1,:));
        end
    end
end
% a = permute(decode_symbol,[2 1 3]);
array=reshape(decode_symbol,M*N,[]);

%-----   scatter plots  ------
lane=1;
figure;
subplot(2,4,1);
scatplot(array(1,:),array(1+M,:))
subplot(2,4,2);
scatplot(array(2,:),array(2+M,:))
subplot(2,4,3);
scatplot(array(3,:),array(3+M,:))

subplot(2,4,4);
scatplot(array(4,:),array(4+M,:))

subplot(2,4,5);
scatplot(array(5,:),array(5+M,:))

subplot(2,4,6);
scatplot(array(6,:),array(6+M,:))

subplot(2,4,7);
scatplot(array(7,:),array(7+M,:))

subplot(2,4,8);
scatplot(array(8,:),array(8+M,:))
% heatscatterplot(array(lane,:)+sqrt(-1)*array(lane+M,:));
%-----   BER and GMI  ------
for c=1:1:M*N
    ber(c)=theore_ber(symbol(c:M*N:end),array(c:M*N:end),1);
end
BER_total=sum(ber)/(M*N);

BER_total=BER_total'

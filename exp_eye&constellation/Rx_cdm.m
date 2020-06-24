clear;clc;close all;
symbolrate=54;
% symbolrate=65;
samplerate=160;
M = 8;%%%%%%%%%% CDM lane
% M = 16;%%%%%%%%%% CDM lane
Walshcode=hadamard(M);
%---- this file for data receiving and processing for experimental use
load('cdm8.mat');
% load('cdm16.mat');

load cdm8_54(a)0.mat wav
% load cdm8_54(b)0.mat wav
% load cdm8_54_BI_(a)0.mat wav
% load cdm8_54_BI_(b)0.mat wav

% save rx_dso16_65(a).mat wav
% save rx_dso16_65(b).mat wav
% save rx_dso16_65_BI(a).mat wav
% save rx_dso16_65_BI(b).mat wav

%-------- processing use  ------------
%------------------------
point=start(symbolrate, samplerate, Tx_data.^2, wav);

wav=wav(point:length(wav)/5);
wav = downsample_data(wav, symbolrate, samplerate, 0.1);
% figure;plot(wav);hold on;plot(Tx_data)

wav = LMS_FFE(Tx_data.^2, wav, 121, 0.0001);
% figure;plot(wav);hold on;plot(Tx_data.^2)

%-----   Rx processing for each frame ------
Rx_symbol = wav;
Rx_symbol=Rx_symbol(1:floor(length(Rx_symbol)/M)*M);
for c=1:1:M
    decode_symbol(c,:)=ds_demod(Walshcode(c,:)',Rx_symbol);
end
decode_symbol=reshape(decode_symbol,M,[]);
%-----   scatter plots  ------

for lane=1:1:M
eyenb(decode_symbol(lane,1:20000),160);
end
%-----   BER and GMI  ------
L=13;
binary=idinput(2^L -1, 'prbs'); 
symbol=Graymap(binary,2);

for c=1:1:M
    ber(c)=theore_ber(symbol(c:M:end),decode_symbol(c:M:end),1);
end
BER_total=sum(ber)/(M);

BER_total=BER_total'

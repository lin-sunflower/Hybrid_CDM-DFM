clear;clc;close all;
upsampling=5;
band=2;
L=13;
m=16;
Data_B=(idinput(2^L -1, 'prbs')+1)/2; 
Data_S=qammod(Data_B(1:floor(length(Data_B)/log2(m))*log2(m)),m,'InputType','bit','UnitAveragePower',true);
% scatterplot(Data_S)
a=0.1;                %% Upsampling
Fd=1/upsampling;
Fs=1;
L=20;
w=rcosine(Fd, Fs, 'fir/sqrt', a, L);
x=-L:Fd:L;
f1_i=w.*sin((2*1-1+0.1)*pi*x);
f1_q=w.*cos((2*1-1+0.1)*pi*x);
f2_i=w.*sin((2*2-1+0.25)*pi*x);
f2_q=w.*cos((2*2-1+0.25)*pi*x);
%-------------   orthogonality plots ------------
% c12=conv(f1,f2);
% c11=conv(f1,f1);
% c22=conv(f2,f2);
% figure;plot(c11);hold on;plot(c22);hold on;plot(c12);hold on;plot(zeros(1,length(c12)))
% % figure;plot(f1);hold on;plot(f2)
% figure;plot(10*log10(abs(fft(c11))));hold on;plot(10*log10(abs(fft(c22))));
% hold on;plot(10*log10(abs(fft(conv(w,w)))))
%-------------   
k=zeros(1,Fs/Fd-1);k=[1 k]';%%上采样函数
Data_S1=Data_S(1:2:end-1);
Data_S2=Data_S(2:2:end);
up_idata1=kron(real(Data_S1),k);
up_qdata1=kron(imag(Data_S1),k);
up_idata2=kron(real(Data_S2),k);
up_qdata2=kron(imag(Data_S2),k);
%-------------   
cap=filter(f1_i,1,up_idata1)+filter(f1_q,1,up_qdata1)+filter(f2_i,1,up_idata2)+filter(f2_q,1,up_qdata2);
figure;plot(10*log10(abs(fft(cap))))

de_idata1=filter(f1_i,1,cap);
de_qdata1=filter(f1_q,1,cap);
de_idata2=filter(f2_i,1,cap);
de_qdata2=filter(f2_q,1,cap);

dsample_idata1=de_idata1(length(w):Fs/Fd:end);
dsample_qdata1=de_qdata1(length(w):Fs/Fd:end);
dsample_idata2=de_idata2(length(w):Fs/Fd:end);
dsample_qdata2=de_qdata2(length(w):Fs/Fd:end);

scatterplot(dsample_idata1+sqrt(-1)*dsample_qdata1)
scatterplot(dsample_idata2+sqrt(-1)*dsample_qdata2)
% figure;plot(dsample_qdata);hold on;plot(imag(Data_S))
% figure;plot(w)
% figure;plot(c11)
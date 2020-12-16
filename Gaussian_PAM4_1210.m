clear;clc;
% close all;
BW_indx =1;
for BW= 25e9:1e9:50e9
    %--------------------------------
    for Length=30
        %     for alpha=0:0.1:0.2
        alpha=0;
        %--------------------------------
        upsampling =1;
        Fd=56e9;     %%% symbol rate in the physical layer
        Fs=Fd*24/upsampling;     %% sampling rate
        H=gaussfir(1/(Fs)*BW,50,Fs/Fd);
        %-------      DA       ------
        a=0.5;
        DA_length=5;
        w_DA=rcosine(Fd, Fs, 'fir/sqrt', a, DA_length);
        k=zeros(1,Fs/Fd-1);k=[1 k];
        % figure;plot(w_DA)
        %--------------------------------
        lw=5e6;     %% linewidth
        c_light=3e8;
        n_1=1.45;
        Td=round(Length*n_1/c_light*Fs);
        wc=2*pi*192e12;
        
        %--------------------------------
        for j=1:1:25
            SNR(j)=j-1;
            snr=10^(SNR(j)/10);
            L=19;
            Binary=idinput(2^L -1, 'prbs');
            Symbol=Graymap(Binary,2);
            %%%%%------- normalized power
            Tx_data = Symbol;
            Tx_data=Tx_data-mean(Tx_data);
            power=sum(Tx_data.^2)/(length(Tx_data));
            Tx_data=Tx_data/sqrt(power);
            %-------      DA       ------
            up_Symbol=kron(Tx_data,k);
            Tx_data=filter(w_DA,1,up_Symbol);
            %--------------------------------
            Tx_data = Tx_data -min(Tx_data);
            Tx_data=sqrt(Tx_data);     %% For square-law detection (MZ)
            N_inten=sqrt(1/snr).*randn(1,length(Tx_data));
            N_phase=sqrt(2*pi*lw*Td/Fs).*randn(length(Tx_data),1);
            %-----------SNR verification------------
            %     sum(N_inten.^2)/(length(N_inten))
            %     Rx_data =Tx_data.^2;
            %     Rx_data = Rx_data -mean(Rx_data);
            %     sum(Rx_data.^2)/(length(Rx_data))
            %%%%%------- BER----------
            for t=1:1:Td
                Rx_data(t)=0;
            end
            for t=Td+1:1:length(Tx_data)
                Rx_data(t)=(1-alpha^2)*Tx_data(t)^2+alpha^2*Tx_data(t-Td)^2 +2*alpha*Tx_data(t)*Tx_data(t-Td)*cos(wc*Td)+N_inten(t);
                %         Rx_data(t)=(1-alpha^2)*Tx_data(t)^2+alpha^2*Tx_data(t-Td)^2 +2*alpha*Tx_data(t)*Tx_data(t-Td)*cos(N_phase(t))+N_inten(t);
            end
            
            %      --------------  bandwidth limitation  ----------------
            K=zeros(1,Fs/Fd-1);K=[1 K];
            Rx_data=kron(Rx_data,K);
            Rx_data=conv(Rx_data,H);
            Rx_data=Rx_data(length(H)/2:Fs/Fd:end-length(H)/2);  %%% 20 is the sample resolution of channel
            %-------      AD       ------
            Rx_data = Rx_data - mean(Rx_data);
            Rx_data=filter(w_DA,1,Rx_data);
            Rx_symbol=Rx_data(length(w_DA):Fs/Fd:end);
            %     figure;plot(Symbol);hold on;plot(Rx_data)
            Start=ceil(Td/(Fs/Fd));
            BER(j)=theore_ber(Symbol,Rx_symbol,Start);
            % GMI(j)=gmi_1d(Symbol,Rx_symbol,4,Start)/2;
        end
        power_penalty(BW_indx)=SNR_penalty(SNR,BER, 3.8e-3)
%         hold on;plot(SNR,BER,'-*')
%         grid on;
    end
    BW_indx = BW_indx+1;
end



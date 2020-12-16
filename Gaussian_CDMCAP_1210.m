clear;clc;
% close all;
M=4;%%%%%%%%%% CDM lane
N=2;%%%%%%%%%% DFM lane
upsampling=3;
%--------------------------------
BW_indx =1;
%--------------------------------
for BW= 25e9:1e9:50e9
    for Length=40
        %     for alpha=0:0.1:0.2
        alpha=0;
        %--------------------------------
        %%%%%------- DFM shaping ----------
        a_cap=0.1;                  %% Upsampling
        cap_length=10;
        Fd_cap=1/3;
        Fs_cap=1;
        w_cap=rcosine(Fd_cap, Fs_cap, 'fir/sqrt', a_cap, cap_length);
        x_cap=-cap_length:Fd_cap/Fs_cap:cap_length;
        f1=w_cap.*sin(1.2*pi*x_cap);
        f2=w_cap.*cos(1.2*pi*x_cap);
        f=[f1;f2];
        k_cap=zeros(1,Fs_cap/Fd_cap-1);k_cap=[1 k_cap]';
        Walshcode=hadamard(M);
        %--------------------------------
        Fd=56e9;     %%% symbol rate in the physical layer
        Fs=Fd*30/upsampling;     %% sampling rate
        H=gaussfir(1/(Fs)*BW,50,Fs/Fd);
        %-------      DA       ------
        a=0.5;
        DA_length=10;
        w_DA=rcosine(Fd, Fs, 'fir/sqrt', a, DA_length);
        k=zeros(1,Fs/Fd-1);k=[1 k];
        % figure;plot(w_DA)
        %--------------------------------
        lw=5e6;     %% 5linewidth
        c_light=3e8;
        n_1=1.45;
        Td=round(Length*n_1/c_light*Fs);
        wc=2*pi*192e12;
        %--------------------------------
        for j=1:1:20
            SNR(j)=j-1;
            snr=10^(SNR(j)/10);
            L=21;
            Binary=idinput(2^L -1, 'prbs');
            symbol=Graymap(Binary,2);
            %-------   CDM coding  ------
            for c1=1:1:N
                for c2=1:1:M
                    code_symbol_I(c2,:)=ds_mod(Walshcode(c2,:)',symbol((c1-1)*M+c2:M*N:end));
                end
                Total_code_I=sum(code_symbol_I);
                up_Symbol=kron(Total_code_I,k_cap');
                filter_data(c1,:)=(-1)^c1*filter(f(c1,:),1,up_Symbol);
            end
            Tx_data=sum(filter_data);
            %%%%%------- normalized power
            Tx_data=Tx_data-mean(Tx_data);
            power=sum(Tx_data.^2)/(length(Tx_data));
            Tx_data=Tx_data/sqrt(power);
            %-------      DA       ------
            up_Symbol=kron(Tx_data,k);
            Tx_data=filter(w_DA,1,up_Symbol);
            %--------------------------------
            Tx_data = Tx_data -min(Tx_data);
            Tx_data=sqrt(Tx_data);     %% For square-law detection (MZ)
            N_phase=sqrt(2*pi*lw*Td/Fs).*randn(length(Tx_data),1);
            N_inten=sqrt(1/snr).*randn(1,length(Tx_data));
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
                %                 Rx_data(t)=(1-alpha^2)*Tx_data(t)^2+alpha^2*Tx_data(t-Td)^2 +2*alpha*Tx_data(t)*Tx_data(t-Td)*cos(N_phase(t))+N_inten(t);
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
            %--------------------------------
            Rx_symbol=Rx_symbol(1:floor(length(Rx_symbol)/M)*M);
            Start=ceil(Td/(Fs/Fd));
            %--------------------------------
            for c1=1:1:N
                defilter_data=filter(f(c1,:),1,Rx_symbol);
                defilter_S=defilter_data(length(w_cap):Fs_cap/Fd_cap:end);
                defilter_symbol(c1,:)=defilter_S(1:floor(length(defilter_S)/M)*M);
                for c2=1:1:M
                    decode_symbol((c1-1)*M+c2,:)= ds_demod(Walshcode(c2,:)',defilter_symbol(c1,:));
                end
            end
            
            for c=1:1:M*N
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %     figure;plot(symbol(c:M*N:end));hold on;plot(decode_symbol(c:M*N:end))
                ber(c)=theore_ber(symbol(c:M*N:end),decode_symbol(c:M*N:end),Start);
                %         gmi(c)=gmi_1d(symbol(c:M*N:end),decode_symbol(c:M*N:end),4,Start);
            end
            BER_total(j)=sum(ber)/(M*N);
            % GMI(j)=gmi_1d(Symbol,Rx_symbol,4,Start)/2;
        end
        %         hold on;plot(SNR,BER_total,'-*')
        %         grid on;
        power_penalty(BW_indx)=SNR_penalty(SNR,BER_total, 3.8e-3)
    end
    BW_indx = BW_indx+1;
end

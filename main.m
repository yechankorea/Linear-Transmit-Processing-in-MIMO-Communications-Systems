clear ;
clc;
close;

%% LINEAR TRANSMIT PROCESSING IN MIMO COMMUNICATIONS SYSTEMS
tic;
example = "fig3" 
    
%basis parameter
    
%% fig2 Receive and transmit filters : MSE e versus SNR r for spatially white noise
if example == "fig2"
    SP.Nt= 2 ;
    SP.Nr = 2;
    SP.H_type  = "Rayleigh" ;


    B = 2;
    R_s = eye(B);
    E_tr = 2;

    experiments = 100000 ;
    P_tx_matrix = eye(SP.Nr);
    G_rx_matrix = eye(SP.Nt);

    SNR_dB = -10:1:30 ;
    SNR = 10.^((SNR_dB./10));

    MSE_RxMF =zeros(1,length(SNR_dB));
    MSE_RxZF =zeros(1,length(SNR_dB));
    MSE_RxWF =zeros(1,length(SNR_dB));
    MSE_TxMF =zeros(1,length(SNR_dB));
    MSE_TxZF =zeros(1,length(SNR_dB));
    MSE_TxWF =zeros(1,length(SNR_dB));

    % SNR = ( E_tr / B ) / (trace(R_n) / Nr) 
    % trace(R_n)= sigma_square * Nr
    sigma_square = E_tr ./ (SNR.*B) ;
   
        for i= 1: length(SNR_dB)
            for expr = 1: experiments 
            R_n= sigma_square(i) .* eye(SP.Nr);
            
            [G_MF,MSE_RxMF_tmp,SNR_receive] = RxMF(SP,P_tx_matrix,R_s,R_n);         
            MSE_RxMF(i) = MSE_RxMF_tmp + MSE_RxMF(i) ;              %RxMF
            
            [G_ZF,MSE_RxZF_tmp,SNR_receive] = RxZF(SP,P_tx_matrix,R_s,R_n);
            MSE_RxZF(i) = MSE_RxZF_tmp + MSE_RxZF(i) ;              %RxZF
            
            [G_WF,MSE_RxWF_tmp,SNR_receive] = RxWF(SP,P_tx_matrix,R_s,R_n,B);
            MSE_RxWF(i) = MSE_RxWF_tmp + MSE_RxWF(i) ;              %RxWF
           
            [P_MF,MSE_TxMF_tmp,SNR_receive] = TxMF(SP,G_rx_matrix,R_s,R_n,E_tr);
            MSE_TxMF(i) = MSE_TxMF_tmp + MSE_TxMF(i) ;              %TxMF
            
            [P_ZF,MSE_TxZF_tmp,SNR_receive] = TxZF(SP,G_rx_matrix,R_s,R_n,E_tr);
            MSE_TxZF(i) = MSE_TxZF_tmp + MSE_TxZF(i) ;              %TxZF
            
            [P_WF,MSE_TxWF_tmp,SNR_receive] = TxWF(SP,G_rx_matrix,R_s,R_n,B,E_tr);
            MSE_TxWF(i) = MSE_TxWF_tmp + MSE_TxWF(i) ;              %TxWF



            end
        fprintf('SNR index: %d \t Elapsed: %.1f s (%.1f hours) \n',i,toc,(toc/3600))
   end

MSE_RxMF = MSE_RxMF./experiments ;
MSE_RxZF = MSE_RxZF./experiments ;
MSE_RxWF = MSE_RxWF./experiments ;
MSE_TxMF = MSE_TxMF./experiments ;
MSE_TxZF = MSE_TxZF./experiments ;
MSE_TxWF = MSE_TxWF./experiments ;

  plot(SNR_dB,MSE_RxMF,'-r')
  hold on
  plot(SNR_dB,MSE_RxZF,'-r+')
  plot(SNR_dB,MSE_RxWF,'-ro')
  
  plot(SNR_dB,MSE_TxMF,'-b')
  plot(SNR_dB,MSE_TxZF,'-b+')
  plot(SNR_dB,MSE_TxWF,'-bo')
  
  hold off
  grid on
  grid minor
  legend('RxMF','RxZF','RxWF','TxMF','TxZF','TxWF')
  xlabel('SNR (in dB)') 
  ylabel('MSE') 
  xlim([-10 30]) 
 set(gca, 'YScale', 'log')
end

%% fig3 : Receive and transmit filters: BER versus SNR for spatially white noise
if example == "fig3"


    SP.Nt= 2 ;
    SP.Nr = 2;
    SP.H_type  = "Rayleigh" ;


    B = 2;
    R_s = eye(B);
    E_tr = 2;

    experiments = 100000 ; % 100,000
    P_tx_matrix = eye(SP.Nr);
    G_rx_matrix = eye(SP.Nt);

    SNR_dB = -10:1:30 ;
    SNR = 10.^((SNR_dB./10));
    sigma_square = E_tr ./ (SNR.*B) ;

    BER_RxMF =zeros(1,length(SNR_dB));
    BER_RxZF =zeros(1,length(SNR_dB));
    BER_RxWF =zeros(1,length(SNR_dB));
    BER_TxMF =zeros(1,length(SNR_dB));
    BER_TxZF =zeros(1,length(SNR_dB));
    BER_TxWF =zeros(1,length(SNR_dB));
    BER_matrix = zeros(length(SNR_dB),6);

    % SNR = ( E_tr / B ) / (trace(R_n) / M) 
    % trace(G * R_n * G')= sigma_square * Nr

    % SNR = ( E_tr / B ) / (trace(G * R_n * G') / M) 
    % trace(G * R_n * G')= sigma_square * trace(G * G')

  
   
        for i= 1: length(SNR_dB)
            for expr = 1: experiments 
            R_n= sigma_square(i) .* eye(SP.Nr);

            % ber_count(SP,R_s,R_n,E_tr,i)
          %  BER_RxMF(i) =BER_RxMF(i)+

%             [G_MF,MSE_RxMF_tmp,SNR_receive] = RxMF(SP,P_tx_matrix,R_s,R_n);    
            BER_matrix(i,:) = ber_count(SP,R_s,R_n,E_tr,i)+ BER_matrix(i,:);
%             BER_RxMF(i) = ber_count('mode_RxMF',SP,R_s,R_n,E_tr)+BER_RxMF(i);
%             BER_RxZF(i) = ber_count('mode_RxZF',SP,R_s,R_n,E_tr)+BER_RxZF(i);
%             BER_RxWF(i) = ber_count('mode_RxWF',SP,R_s,R_n,E_tr)+BER_RxWF(i);
%             BER_TxMF(i) = ber_count('mode_TxMF',SP,R_s,R_n,E_tr)+BER_TxMF(i);
%             BER_TxZF(i) = ber_count('mode_TxZF',SP,R_s,R_n,E_tr)+BER_TxZF(i);
%             BER_TxWF(i) = ber_count('mode_TxWF',SP,R_s,R_n,E_tr)+BER_TxWF(i);

%             BER_RxMF(i) = BER_cal(SNR_receive,M,k) + BER_RxMF(i) ;               %RxMF
            
%             [G_ZF,MSE_RxZF_tmp,SNR_receive] = RxZF(SP,P_tx_matrix,R_s,R_n);
%             BER_RxZF(i) = BER_cal(SNR_receive,M,k) + BER_RxZF(i) ;              %RxZF
%             
%             [G_WF,MSE_RxWF_tmp,SNR_receive] = RxWF(SP,P_tx_matrix,R_s,R_n,B);
%             BER_RxWF(i) = BER_cal(SNR_receive,M,k) + BER_RxWF(i) ;              %RxWF
%            
%             [P_MF,MSE_TxMF_tmp,SNR_receive] = TxMF(SP,G_rx_matrix,R_s,R_n,E_tr);
%             BER_TxMF(i) = BER_cal(SNR_receive,M,k) + BER_TxMF(i) ;              %TxMF
%             
%             [P_ZF,MSE_TxZF_tmp,SNR_receive] = TxZF(SP,G_rx_matrix,R_s,R_n,E_tr);
%             BER_TxZF(i) = BER_cal(SNR_receive,M,k) + BER_TxZF(i) ;              %TxZF
%             
%             [P_WF,MSE_TxWF_tmp,SNR_receive] = TxWF(SP,G_rx_matrix,R_s,R_n,B,E_tr);
%             BER_TxWF(i) = BER_cal(SNR_receive,M,k) + BER_TxWF(i) ;              %TxWF

    % SNR = ( E_tr / B ) / (trace(R_n) / M) 
    % trace(R_n)= sigma_square * Nr

            end
        fprintf('SNR index: %d \t Elapsed: %.1f s (%.1f hours) \n',i,toc,(toc/3600))
   end
BER_matrix =BER_matrix./experiments; 
% BER_RxMF = BER_RxMF./experiments ;
% BER_RxZF = BER_RxZF./experiments ;
% BER_RxWF = BER_RxWF./experiments ;
% BER_TxMF = BER_TxMF./experiments ;
% BER_TxZF = BER_TxZF./experiments ;
% BER_TxWF = BER_TxWF./experiments ;


% P_m = (1-1/sqrt(M))*erfc(sqrt(3*k*SNR(i)/ / (2*(M-1)) )); % 100QAM , M=100, k=7
% P_s =1-(1-P_m)^2;
% P_b=(1/k)*P_s;
    % SNR = ( E_tr / B ) / (trace(R_n) / M) 
    % trace(R_n)= sigma_square * Nr


  semilogy(SNR_dB,BER_matrix(:,1),'-r')
  hold on
  semilogy(SNR_dB,BER_matrix(:,2),'-r+')
  semilogy(SNR_dB,BER_matrix(:,3),'-ro')
  semilogy(SNR_dB,BER_matrix(:,4),'-b')
  semilogy(SNR_dB,BER_matrix(:,5),'-b+')
  semilogy(SNR_dB,BER_matrix(:,6),'-bo')

%   semilogy(SNR_dB,BER_RxMF,'-r')
%   hold on
%   semilogy(SNR_dB,BER_RxZF,'-r+')
%   semilogy(SNR_dB,BER_RxWF,'-ro')
%   semilogy(SNR_dB,BER_TxMF,'-b')
%   semilogy(SNR_dB,BER_TxZF,'-b+')
%   semilogy(SNR_dB,BER_TxWF,'-bo')
  


  hold off
  grid on
  grid minor
  legend('RxMF','RxZF','RxWF','TxMF','TxZF','TxWF')
  xlabel('SNR (in dB)') 
  ylabel('BER') 
 
%      set(gca, 'YScale', 'log')
  %   ylim([10^(-3) 0]) 


end





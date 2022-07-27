function [BER_vector] = ber_count(SP,R_s,R_n,E_tr,i)
%%basis parameter
[B,B] = size(R_s);
P = eye(B) ;
G = eye(B) ;

H = Channel(SP);

BER_vector = zeros(1,6);
%% RXMF
% if mode == "mode_RxMF"
P = eye(B) ;
G = eye(B) ;
alpha_1 = 1;
G_MF = alpha_1 .* R_s * P' * H' *R_n^(-1);      %(7)
G = G_MF;
BER_vector(1) = CAL(SP,R_n,G,H,P);
% end
%% RXZF
% if mode == "mode_RxZF"
P = eye(B) ;
G = eye(B) ;
G_ZF = (P'*H'*R_n^(-1)*H*P)^(-1)*P'*H'*R_n^(-1); %(10)
G= G_ZF;
BER_vector(2) = CAL(SP,R_n,G,H,P);
% end
%% RXWF
% if mode == "mode_RxWF"
P = eye(B) ;
G = eye(B) ;
G_WF = ( P'*H'*R_n^(-1)*H*P + R_s^(-1) )^(-1)*P'*H'*R_n^(-1); %(10)
G=G_WF;
BER_vector(3) = CAL(SP,R_n,G,H,P);
% end
%% TXMF
% if mode == "mode_TxMF"
P = eye(B) ;
G = eye(B) ;
beta_TxMF = sqrt(E_tr/ trace(H' * G' * R_s * G * H));
P_MF = beta_TxMF.* H' * G';   %(19)
P = P_MF;
BER_vector(4) = CAL(SP,R_n,G,H,P);
% end
%% TXZF
% if mode == "mode_TxZF"
P = eye(B) ;
G = eye(B) ;
beta_TxZF = sqrt(E_tr/ ( trace( (G * H * H' * G')^(-1) * R_s )  ) ); %(26)
P_ZF = beta_TxZF .* H' * G' * (G * H * H' * G')^(-1) ; %(25)
P=P_ZF;
BER_vector(5) = CAL(SP,R_n,G,H,P);
% end
%% TXWF
% if mode == "mode_TxWF"
P = eye(B) ;
G = eye(B) ;
F= (H' * G' * G * H ) + ( trace(G * R_n * G')/E_tr ) .*eye(SP.Nt) ;
beta_TxWF = sqrt(E_tr/ (trace(F^(-2) * H' * G' * R_s * G * H) ) ) ; %(39)
P_WF = beta_TxWF .* F^(-1) *H' * G' ; 
P=P_WF;
BER_vector(6) = CAL(SP,R_n,G,H,P);
% end
%%

% %% BER CALCULATOR SETTING
% function [BER] = CAL(SP,R_n,R_s,G,H,P)
% sigma_square = trace(R_n) / SP.Nr ;
% n= (1/sqrt(2))*(sqrt(sigma_square).*rand(2,1)+sqrt(sigma_square).*rand(2,1).*j) ;
% 
% a = rand(1);
% b = rand(1);
% c = rand(1);
% d = rand(1);
% count = 0;
% if a > 0.5
%     a = 1;
% else
%     a = -1;
% end
% if b > 0.5
%     b = 1;
% else
%     b = -1;
% end
% if c > 0.5
%     c = 1;
% else
%     c = -1;
% end
% if d > 0.5
%     d = 1;
% else
%     d = -1;
% end
%     
%     
% s_tx_1 = a + b * j;
% s_tx_2 = c + d * j;
% s_tx = transpose([s_tx_1 s_tx_2]);
% 
% 
% s_rx = G * H * P * s_tx + G * n;
% s_rx = transpose(s_rx);
% s_rx_1 = s_rx(1);
% s_rx_2 = s_rx(2);
% 
% if (( (a<0) && (real(s_rx_1)) >= 0) || ( (a>=0) && (real(s_rx_1) < 0) ) )
%     count=count+1;
% end
% if (( (b<0) && (imag(s_rx_1)) >= 0) || ( (b>=0) && (imag(s_rx_1) < 0) ) )
%     count=count+1;
% end
% if (( (c<0) && (real(s_rx_2)) >= 0) || ( (c>=0) && (real(s_rx_2) < 0) ) )
%     count=count+1;
% end
% if (( (d<0) && (imag(s_rx_2)) >= 0) || ( (d>=0) && (imag(s_rx_2) < 0) ) )
%     count=count+1;
% end
% 
% BER= count / 4;
% 
% end
end

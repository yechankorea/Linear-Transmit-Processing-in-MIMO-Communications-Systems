function [G_ZF,MSE_RxZF,SNR_receive] = RxZF(SP, P, R_s, R_n)

H = Channel(SP);

G_ZF = (P'*H'*R_n^(-1)*H*P)^(-1)*P'*H'*R_n^(-1); %(10)
G= G_ZF;
J_rx = R_s'^(1/2) * P' * H' * R_n^(-1) * H * P *R_s^(1/2); 
 
MSE_RxZF =trace(R_s)- ((trace(R_s))^2 ) / trace( R_s + J_rx^(-1)*R_s); % (16)

SNR_receive = trace(G *H *P * R_s * P' * H' * G') / trace(G * R_n * G');

end

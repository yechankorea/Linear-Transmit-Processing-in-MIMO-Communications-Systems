function [G_MF,MSE_RxMF,SNR_receive] = RxMF(SP, P, R_s, R_n)

H = Channel(SP);
alpha_1 = 1;

G_MF = alpha_1 .* R_s * P' * H' *R_n^(-1);      %(7)
G = G_MF;
J_rx = R_s'^(1/2) * P' * H' * R_n^(-1) * H * P *R_s^(1/2); 
 
MSE_RxMF =trace(R_s)- ((trace(J_rx*R_s))^2 ) / trace( (J_rx^2+J_rx)*R_s); % (15)

SNR_receive = trace(G *H *P * R_s * P' * H' * G') / trace(G * R_n * G');

end

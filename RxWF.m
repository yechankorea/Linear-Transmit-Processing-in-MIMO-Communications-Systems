function [G_WF,MSE_RxWF,SNR_receive] = RxWF(SP, P, R_s, R_n,B)

H = Channel(SP);

G_WF = ( P'*H'*R_n^(-1)*H*P + R_s^(-1) )^(-1)*P'*H'*R_n^(-1); %(10)
G=G_WF;
J_rx = R_s'^(1/2) * P' * H' * R_n^(-1) * H * P *R_s^(1/2); 
 
MSE_RxWF =trace(R_s)- trace( (J_rx+ eye(B))^(-1)*J_rx * R_s ); % (17)

SNR_receive = trace(G *H *P * R_s * P' * H' * G') / trace(G * R_n * G');


end

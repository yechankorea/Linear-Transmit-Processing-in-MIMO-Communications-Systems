function [P_MF,MSE_TxMF,SNR_receive] = TxMF(SP, G, R_s, R_n,E_tr)

H = Channel(SP);

beta_TxMF = sqrt(E_tr/ trace(H' * G' * R_s * G * H));
P_MF = beta_TxMF.* H' * G';   %(19)
P=P_MF;

J_tx = (E_tr/trace(G*R_n*G')) .* (G * H * H' * G') ;

MSE_TxMF =trace(R_s)- ((trace(J_tx*R_s))^2 ) / trace( (J_tx^2+J_tx)*R_s); % (20)
SNR_receive = trace(G *H *P * R_s * P' * H' * G') / trace(G * R_n * G');
end

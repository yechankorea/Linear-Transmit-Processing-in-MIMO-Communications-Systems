function [P_WF,MSE_TxWF,SNR_receive] = TxWF(SP, G, R_s, R_n,B,E_tr)

H = Channel(SP);

F= (H' * G' * G * H ) + ( trace(G * R_n * G')/E_tr ) .*eye(SP.Nt) ;
beta_TxWF = sqrt(E_tr/ (trace(F^(-2) * H' * G' * R_s * G * H) ) ) ; %(39)

P_WF = beta_TxWF .* F^(-1) *H' * G' ; 
P=P_WF;
J_tx = (E_tr/trace(G*R_n*G')) .* (G * H * H' * G') ;
 
MSE_TxWF =trace(R_s)- trace( (J_tx+ eye(B))^(-1)*J_tx * R_s ); % (17)

SNR_receive = trace(G *H *P * R_s * P' * H' * G') / trace(G * R_n * G');
end
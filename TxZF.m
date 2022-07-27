function [P_ZF,MSE_TxZF,SNR_receive] = TxZF(SP, G, R_s, R_n,E_tr)

H = Channel(SP);
beta_TxZF = sqrt(E_tr/ ( trace( (G * H * H' * G')^(-1) * R_s )  ) ); %(26)

P_ZF = beta_TxZF .* H' * G' * (G * H * H' * G')^(-1) ; %(25)
P=P_ZF;
J_tx = (E_tr/trace(G*R_n*G')) .* (G * H * H' * G') ;
 
MSE_TxZF =trace(R_s)- ((trace(R_s))^2 ) / trace( R_s + J_tx^(-1)*R_s); % (30)

SNR_receive = trace(G *H *P * R_s * P' * H' * G') / trace(G * R_n * G');

end

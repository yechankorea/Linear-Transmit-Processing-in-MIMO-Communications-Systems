function [BER] = CAL(SP,R_n,G,H,P)
sigma_square = trace(R_n) / SP.Nr ;
n= (1/sqrt(2)).*(sqrt(sigma_square).*randn(2,1)+sqrt(sigma_square).*randn(2,1).*j) ;
%% symbol generator B=2  
count = 0;
BER = 0;
for time_slot = 1: 100
a=rand(1);
b=rand(1);
c=rand(1);
d=rand(1);
if  a>= 0.5
    a = 1;
else
    a=-1;
end
if  b>= 0.5
    b = 1;
else
    b=-1;
end
if  c>= 0.5
    c = 1;
else
    c=-1;
end
if  d>= 0.5
    d = 1;
else
    d=-1;
end

% B_1 = rand(1);
% B_2 = rand(1);
% 
% if B_1 >= 0.75
%     bit_symbol_1 = '11';
% elseif B_1>= 0.5
%     bit_symbol_1 = '10';
% elseif B_1>= 0.25
%     bit_symbol_1 = '01';
% else
%     bit_symbol_1 = '00';
% end
% 
% if B_2 >= 0.75
%     bit_symbol_2 = '11';
% elseif B_2>= 0.5
%     bit_symbol_2 = '10';
% elseif B_2>= 0.25
%     bit_symbol_2 = '01';
% else
%     bit_symbol_2 = '00';
% end
%     %% symbol encoding
%    if bit_symbol_1 == '11'
%        symbol_map_1 = 1+j ;
%    elseif bit_symbol_1 == '10'
%        symbol_map_1 = 1-j ;
%    elseif bit_symbol_1 == '01'
%        symbol_map_1 = -1+j ;
%    elseif bit_symbol_1 =='00'
%        symbol_map_1 = -1-j ;
%    end
% 
%    if bit_symbol_2 == '11'
%        symbol_map_2 = 1+j ;
%    elseif bit_symbol_2 == '10'
%        symbol_map_2 = 1-j ;
%    elseif bit_symbol_2 == '01'
%        symbol_map_2 = -1+j ;
%    elseif bit_symbol_2 =='00'
%        symbol_map_2 = -1-j ;
%    end
%     symbol_map_tx = transpose( [symbol_map_1, symbol_map_2]);
%     %% channel
%     symbol_map_rx = G * H * P * symbol_map_tx + G * n;  
%     %% decoding
%     symbol_map_rx = transpose(symbol_map_rx);
%     
%     if ( (real(symbol_map_rx(1)) >= 0) && (imag(symbol_map_rx(1))>=0) )
%         bit_symbol_rx_1 = '11';
%     
%     elseif ( (real(symbol_map_rx(1)) >= 0) && (imag(symbol_map_rx(1))<0)        
%         bit_symbol_rx_1 = '10';
%     
%     elseif ( (real(symbol_map_rx(1)) < 0) && (imag(symbol_map_rx(1))>=0) )
%         bit_symbol_rx_1 = '01';
% 
%     elseif ( (real(symbol_map_rx(1)) < 0) && (imag(symbol_map_rx(1))<0)        
%         bit_symbol_rx_1 = '00';
%     end
% 
%      if ( (real(symbol_map_rx(2)) >= 0) && (imag(symbol_map_rx(2))>=0) )
%         bit_symbol_rx_1 = '11';
%     
%      elseif ( (real(symbol_map_rx(2)) >= 0) && (imag(symbol_map_rx(2))<0)        
%         bit_symbol_rx_1 = '10';
%     
%      elseif ( (real(symbol_map_rx(2)) < 0) && (imag(symbol_map_rx(2))>=0) )
%         bit_symbol_rx_1 = '01';
%     
%      elseif ( (real(symbol_map_rx(2)) < 0) && (imag(symbol_map_rx(2))<0)        
%         bit_symbol_rx_1 = '00';
%     end
% 
%     bit_symbol_rx = transpose([bit_symbol_rx_1, bit_symbol_rx_2]);



    
s_tx_1 = a + b * j;
s_tx_2 = c + d * j;
s_tx = transpose([s_tx_1 s_tx_2]);


s_rx = G * H * P * s_tx + G * n;
s_rx = transpose(s_rx);
s_rx_1 = s_rx(1);
s_rx_2 = s_rx(2);

if (( (a<0) && (real(s_rx_1)) >= 0) || ( (a>=0) && (real(s_rx_1) < 0) ) )
    count=count+1;
end
if (( (b<0) && (imag(s_rx_1)) >= 0) || ( (b>=0) && (imag(s_rx_1) < 0) ) )
    count=count+1;
end
if (( (c<0) && (real(s_rx_2)) >= 0) || ( (c>=0) && (real(s_rx_2) < 0) ) )
    count=count+1;
end
if (( (d<0) && (imag(s_rx_2)) >= 0) || ( (d>=0) && (imag(s_rx_2) < 0) ) )
    count=count+1;
end

end
BER= count/400;
end
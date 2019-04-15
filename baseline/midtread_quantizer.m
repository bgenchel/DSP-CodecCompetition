function [ret_value] = midtread_quantizer(x, R)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    MIDTREAD_QUANTIZER     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q = 2 / (2^R - 1);      
q = quant(x,Q);
s = q<0;    
ret_value = uint16(abs(q)./Q + s*2^(R-1));
end
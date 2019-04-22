function [x, bitsleft] = allocate_uni(y, bits_num, N, Fs, num_subbands, numBandsToIgnore)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        ALLOCATE           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x=allocate(y,b,sb,N,Fs)
% Allocates b bits to the 25 subbands
% of y (a length N/2 MDCT, in dB SPL)
%
% modified by Jiawen: allocate same number of bits to positive 'bits'

bits_per_frame = bits_num;
        
bits(floor(bark( (Fs/2)*(1:N/2)/(N/2) )) +1) = 0;

for ii=1:N/2
    bits(floor(bark( (Fs/2)*ii/(N/2) )) +1) = max(bits(floor(bark( (Fs/2)*ii/(N/2) )) +1) , ceil( y(ii)/6 ));
end

indices = find(bits(1:end) < 2);
bits(indices(1:end)) = 0;
bits(end-numBandsToIgnore:end) = 0; % artifact reduction

% NEED TO CALCULATE SAMPLES PER SUBBAND
n = 0:N/2-1;
f_Hz = n*Fs/N;
f_kHz = f_Hz / 1000;
z = 13*atan(0.76*f_kHz) + 3.5*atan((f_kHz/7.5).^2);  % *** bark frequency scale
crit_band = floor(z)+1;
num_crit_bands = max(crit_band);
num_crit_band_samples = zeros(num_crit_bands,1);
for ii=1:N/2
    num_crit_band_samples(crit_band(ii)) = num_crit_band_samples(crit_band(ii)) + 1;
end

%{
x = zeros(1,num_subbands);
bitsleft = bits_per_frame;
[~,ii]=max(bits);
while bitsleft > num_crit_band_samples(ii)
    [~,ii]=max(bits); % this line should be at the end of the while loop...
    x(ii) = x(ii) + 1;
    bits(ii) = bits(ii) - 1;
    bitsleft=bitsleft-num_crit_band_samples(ii);
end
end
%}

% if we allocate the same number of bits for each band...
t = floor(bits_per_frame/sum(num_crit_band_samples(bits>0)));
if t > 15
    t = 15;
end
x = zeros(1, num_subbands);
indices = find(bits > 0);
for iii = 1:length(indices)
    x(indices(iii)) = t;
end
bitsleft =  bits_per_frame - x * num_crit_band_samples;

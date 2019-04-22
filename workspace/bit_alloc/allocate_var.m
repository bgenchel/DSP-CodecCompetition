function [x, bitsleft] = allocate_var(y, bits_num, N, Fs, num_subbands, numBandsToIgnore)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        ALLOCATE           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x=allocate(y,b,sb,N,Fs)
% Allocates b bits to the 25 subbands
% of y (a length N/2 MDCT, in dB SPL)

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
var_subbands = zeros(1, num_subbands);
for ii=1:num_subbands
    indices = find(crit_band==ii);
    energy_band = y(indices);
    if max(energy_band) > 0
        energy_band = energy_band/max(energy_band);
    end
    var_subbands(ii) = std(energy_band);
end
indices = var_subbands==0;

x = zeros(1,num_subbands);

bitsleft = bits_per_frame - sum(num_crit_band_samples(indices));
x(indices) = 1;

var_subbands(indices) = Inf;

[~,ii]=min(var_subbands);
while bitsleft > num_crit_band_samples(ii)
    if(x(ii) < 15)
        x(ii) = x(ii) + 1;
        bitsleft=bitsleft-num_crit_band_samples(ii);
    end
    var_subbands(ii) = var_subbands(ii) + 1;
    [~,ii]=min(var_subbands); % this line should be at the end of the while loop...
end
end
function [x, bitsleft] = allocate_default(y, bits_num, N, Fs, num_subbands, numBandsToIgnore)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     ALLOCATE DEFAULT      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the default allocate method (same as the baseline, but with
% different input.) All bits_num bits are assigned to the data. No bits 
% for 'scale bits'.
%
% Forget about the allocate function in the baseline folder...
% input:    y, dB SPL for one frame
%           bits_num, number of bits assigned for this frame
%           N, FFT size
%           Fs, sample rate

bits_per_frame = bits_num; % floor(((bitrate/Fs)*(N/2)) - (scalebits*num_subbands));

% bits: used as reference when assigning bits
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

x = zeros(1,num_subbands);
bitsleft = bits_per_frame;
[~,ii]=max(bits);
while bitsleft > num_crit_band_samples(ii)
    if(x(ii) < 15)
        x(ii) = x(ii) + 1;
        bitsleft=bitsleft-num_crit_band_samples(ii);
    end
    bits(ii) = bits(ii) - 1;
    [~,ii]=max(bits); % this line should be at the end of the while loop...
end
end

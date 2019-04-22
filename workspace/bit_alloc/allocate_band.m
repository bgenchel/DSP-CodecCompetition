function [x, bitsleft] = allocate_band(y, bits_num, N, Fs, num_subbands, numBandsToIgnore, method)
x = zeros(1, num_subbands);
bitsleft = 0;
% choose the method to allocate bits to bands
if method == 'default' % original method (bug fixed
    x = allocate_default(y, bits_num, N, Fs, num_subbands, numBandsToIgnore);
end
if method == 'uniform'
    x = allocate_uni(y, bits_num, N, Fs, num_subbands, numBandsToIgnore);
end
if method == 'var    '
    x = allocate_var(y, bits_num, N, Fs, num_subbands, numBandsToIgnore);
end
if method == 'altn   '
    [x, bitsleft] = allocate_altn(y, bits_num, N, Fs, num_subbands, numBandsToIgnore);
end
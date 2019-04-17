function [bit_alloc, gain, data] = allocate_encode_all(y, bitrate, scalebits, N, Fs, frames)
% allocate bits
% encode
% return bit_alloc, gain, data for all frames
frame_num = size(y, 1);
num_subbands = floor(fftbark(N/2,N/2,Fs))+1;

% initialize output
bit_alloc = zeros(frame_num, num_subbands);
gain = zeros(num_subbands, frame_num);
data = zeros(N/2, frame_num);

% allocate bits per frame
% default: same number of bits per frame
bit_per_frame = ones(1, frame_num) * bitrate;

for frame_count=1:frame_num
    bit_alloc(frame_count,:) = allocate_uni(y(frame_count,:),bit_per_frame(frame_count),scalebits,N,Fs);
    [gain(:,frame_count),data(:,frame_count)] = p_encode(mdct(frames(frame_count,:)),Fs,N,bit_alloc(frame_count,:),scalebits);
end
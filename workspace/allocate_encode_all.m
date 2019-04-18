function [bit_alloc, gain, data] = allocate_encode_all(y, bitrate, scalebits, N, Fs, frames)
% allocate bits
% encode
% return bit_alloc, gain, data for all frames
frame_num = size(y, 1);

% allocate
bit_alloc = allocate_main(y, bitrate, scalebits, N, Fs, "default_var");

for frame_count=1:frame_num
    % bit_alloc(frame_count,:) = allocate_band(y(frame_count,:),bit_per_frame(frame_count),N,Fs, num_subbands, bandsIgnore, 'default');
    [gain(:,frame_count),data(:,frame_count)] = p_encode(mdct(frames(frame_count,:)),Fs,N,bit_alloc(frame_count,:),scalebits);
end
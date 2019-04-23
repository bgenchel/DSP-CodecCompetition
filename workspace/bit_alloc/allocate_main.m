function bit_alloc = allocate_main(y, bitrate, scalebits, N, Fs, method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      ALLOCATE MAIN        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assign bits to each frame
% bits_frame: (1 x frame_num)

frame_num = size(y, 1);
num_subbands = floor(fftbark(N/2,N/2,Fs))+1;
bandsIgnore = 2*floor(128000/bitrate);
total_bits = floor(((bitrate/Fs)*(N/2) - (scalebits*num_subbands))*frame_num);

% Intuitive ways of allocation

% equal number of bits for every frame:
if method == 'default_default' % default bit allocation
    bits_average = floor(((bitrate/Fs)*(N/2)) - (scalebits*num_subbands));
    bits_frame = ones(1, frame_num) * bits_average;
    bitsleft = 0;
    for frame_count=1:frame_num
        [bit_alloc(frame_count,:), bitsleft] = allocate_band( ... 
            y(frame_count, :), ... 
            bits_frame(frame_count) + bitsleft, ...
            N, Fs, num_subbands, bandsIgnore, 'default' ...
        );
    end
end

if method == 'default_uniform' % uniformly allocate to bands in a frame
    bits_average = floor(((bitrate/Fs)*(N/2)) - (scalebits*num_subbands));
    bits_frame = ones(1, frame_num) * bits_average;
    bitsleft = 0;
    for frame_count=1:frame_num
        [bit_alloc(frame_count, :), bitsleft] = allocate_band( ...
            y(frame_count, :), ... 
            bits_frame(frame_count) + bitsleft, ...
            N,Fs, num_subbands, bandsIgnore, 'uniform' ...
        );
    end
end

if method == 'default_var' % uniformly allocate to bands in a frame
    bits_average = floor(((bitrate/Fs)*(N/2)) - (scalebits * num_subbands));
    bits_frame = ones(1, frame_num) * bits_average;
    bitsleft = 0;
    for frame_count=1:frame_num
        [bit_alloc(frame_count, :), bitsleft] = allocate_band( ...
            y(frame_count, :), ...
            bits_frame(frame_count) + bitsleft, ... 
            N, Fs, num_subbands, bandsIgnore, 'var    ' ...
        );
    end
end

% different number of bits for every frame

% not good
if method == 'energy_default' % assign by log energy across time, default/uniform within frame
    % although logarithmic addition doesn't really make sense...let's try
    sum_y = sum(y, 2);
    p = sum_y/sum(sum_y);
    bits_frame = floor(p * total_bits);
    for frame_count=1:frame_num
        bit_alloc(frame_count,:) = allocate_band(y(frame_count,:),bits_frame(frame_count),N,Fs, num_subbands, bandsIgnore, 'default');
    end
end

if method == 'maxeng_default' 
    % assign by log energy across time, default/uniform within frame
    % although logarithmic addition doesn't really make sense ... let's try
    sum_y = max(y, [], 2);
    p = sum_y/sum(sum_y);
    bits_frame = floor(p * total_bits);
    for frame_count=1:frame_num
        bit_alloc(frame_count,:) = allocate_band(y(frame_count,:),bits_frame(frame_count),N,Fs, num_subbands, bandsIgnore, 'default');
    end
end

if method == 'sqrtstd_default'
    t = sqrt(std(y, 0, 2));
    p = t./sum(t);
    bits_frame = floor(p * total_bits);
    for frame_count=1:frame_num
        bit_alloc(frame_count,:) = allocate_band(y(frame_count,:),bits_frame(frame_count),N,Fs, num_subbands, bandsIgnore, 'default');
    end
end

if method == 'sqrtstd_uniform'
    t = sqrt(std(y, 0, 2));
    p = t./sum(t);
    bits_frame = floor(p * total_bits);
    for frame_count=1:frame_num
        bit_alloc(frame_count,:) = allocate_band(y(frame_count,:),bits_frame(frame_count),N,Fs, num_subbands, bandsIgnore, 'uniform');
    end
end
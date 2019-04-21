function decodedFile = full_codec_bitalloc(originalFile, bitrate, decodedFile, codedFile)
% FULL_CODEC encodes and decodes an audio file
%   decodedFile = FULL_CODEC(originalFile) encodes and decodes original file

%   decodedFile = FULL_CODEC(originalFile,bitrate) lets you specify the
%       encoded bitrate in bits per second (128000 default)

%   FULL_CODEC(originalFile,bitrate,decodedFile) lets you specify the
%       decoded file

%   FULL_CODEC(originalFile,bitrate,decodedFile,codedFile) also lets you
%       specify the encoded file
%
%   This code is based loosely on the assignments in the textbook:
%   Bosi, Marina, and Richard E. Goldberg. "Introduction to digital audio
%     coding and standards". Springer, 2003.
% 
%   See also AUDIOWRITE

% Github location:
%   https://github.com/JonBoley/SimpleCodec.git
% Bug reports:
%   https://github.com/JonBoley/SimpleCodec/issues
% 
% Copyright 2002, Jon Boley (jdb@jboley.com)

if nargin < 4
    codedFile = 'audio/yourfile.jon';
end
if nargin < 3
    decodedFile = 'yourfile_decoded.wav';
end
if nargin < 2
    bitrate = 64000; % Required by the competition
end
if nargin < 1
    originalFile = 'yourfile.wav';
end

scalebits = 4;
N = 2048; % framelength

[Y, Fs] = audioread(originalFile);
Y = Y(:, 1); % just use the first channel

sig = sin(2 * pi * 1000 * (1:N/2) / Fs);
win = (0.5 - 0.5 * cos((2 * pi * (1:N/2) - 0.5) / (N/2))); 
fftmax = max(abs(fft(sig .* win))); % defined as 96dB, strongest frequency

%   Enframe Audio (divide into frames) 
frames = enframe(Y, N, N/2); 

% Write File Header 
fid = fopen(codedFile, 'w');
fwrite(fid, Fs, 'ubit16'); % Sampling Frequency
fwrite(fid, N, 'ubit12'); % Frame Length
fwrite(fid, bitrate, 'ubit18'); % Bit Rate
fwrite(fid, scalebits, 'ubit4'); % Number of Scale Bits per Sub-Band
fwrite(fid, length(frames(:,1)), 'ubit26'); % Number of frames
    
% fftbark - converts fft bin number to bark scale
% bark scale maps frequency ranges in Hz to integer value. It is, in
% itself, a means of grouping bins.
numBands = floor(fftbark(N/2, N/2, Fs)) + 1;

% Compute signal spectrum after masking
New_FFT_all = spectrumMasking_m(frames, N, Fs, fftmax);

% Allocate bits and Encode
[bit_alloc_all, Gain_all, Data_all] = allocate_encode_all(New_FFT_all, bitrate, scalebits, N, Fs, frames);
    
% Write Audio Data to File
for frame_count=1:length(frames(:,1))
    qbits = sprintf('ubit%i', scalebits);
    fwrite(fid, Gain_all(:,frame_count), qbits);
    fwrite(fid, bit_alloc_all(frame_count,:), 'ubit4');
    for ii=1:numBands
        indices = find((floor(fftbark(1:N/2,N/2,Fs))+1)==ii);
        qbits = sprintf('ubit%i', bit_alloc_all(frame_count,ii)); % bits(floor(fftbark(i,framelength/2,48000))+1)
        if bit_alloc_all(frame_count,ii)>0
            fwrite(fid, Data_all(indices(1):indices(end),frame_count) ,qbits);
        end
    end
end
fclose(fid);

% RUN DECODER
disp('Decoding...');
p_decode(codedFile,decodedFile);

disp('Okay, all done!');




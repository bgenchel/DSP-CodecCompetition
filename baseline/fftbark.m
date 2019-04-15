function b = fftbark(bin, N, Fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          FFTBARK          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% b=fftbark(bin,N,Fs)
% Converts fft bin number to bark scale
%
% from wikipedia - The Bark scale is a psychoacoustical scale proposed by 
% Eberhard Zwicker in 1961. It is named after Heinrich Barkhausen who 
% proposed the first subjective measurements of loudness. One definition
% of the term is "...a frequency scale on which equal distances correspond 
% with perceptually equal distances. Above about 500 Hz this scale is more 
% or less equal to a logarithmic frequency axis. Below 500 Hz the Bark 
% scale becomes more and more linear.
% 
% The scale ranges from 1 to 24 and corresponds to the first 24 critical 
% bands of hearing.
%
% N is the fft length
% Fs is the sampling frequency

f = bin * (Fs/2) / N;
b = 13 * atan(0.76 * f/1000) + 3.5 * atan((f/7500).^2); 
end
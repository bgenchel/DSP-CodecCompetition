function New_FFT_all = temporalTNSpectrumMasking(frames, N, Fs, fftmax)
%%%%%%%%%%%%%%%%%%%%%%%%
% SPECTRUM MASKING
% w/ temporal masking & 
% tonal / noise masking
%%%%%%%%%%%%%%%%%%%%%%%%
% frames - audio frames 
% N - frame length
% Fs - sampling frequency
% fftmax - reference for 96dB, the 'loudest' frequency, 1KHz
%%

% apply hann window
w = hann(size(frames,2));
for i=1:size(frames,1)
    frames(i,:) = frames(i,:).*w';
end

% initialize output
masks = zeros(length(frames(:,1)), N/2);
spls = zeros(length(frames(:,1)), N/2);
New_FFT_all = zeros(length(frames(:,1)), N/2);

% process each frame
for frame_count=1:length(frames(:, 1))

    % print processing status
    if mod(frame_count, 100) == 0
        outstring = sprintf('Now Masking Frame %i of %i', frame_count, length(frames(:,1)));
        disp(outstring);
    end
    
    % take fft of windowed frame
    fft_frame = fft(frames(frame_count, :));

    % if frame is perfect silence, nothing to be done
    if fft_frame == zeros(1, N)
        New_FFT_all(frame_count,:) = zeros(1, N/2);
    else    
        len = length(fft_frame);

        % transform fft to spl (sound pressure level) values
        fft_spl = 96 + 20 * log10(abs(fft_frame)/fftmax);
        fft_spl = fft_spl(1:N/2);
        
        % frequency value in Hz for each bin in FFT
        f_kHz = (1:Fs/N:Fs/2)/1000;
        
        %{
        tonal = [];
        noise = [];
        
        % Find Peaks
        for i=1:size(f_kHz,2)
            if f_kHz <= 11
                width = 1;
            else
                width = 2;
            end
            st = max(1, i-width);
            ed = min(N/2, i+width);
            f = 1;
            for ii=st:ed
                if ii == i
                    continue;
                end
                if fft_spl(i)-fft_spl(ii) < 0
                    f = 0;
                    break;
                end
                if abs(ii-i) ~= 1 && fft_spl(i)-fft_spl(ii) <= 7
                    % noise
                    noise = [noise i];
                    f = 0;
                    break;
                end
            end
            if f == 1
                tonal = [tonal i];
            end
        end
        %}
        
        % Find Peaks
        centers = find(diff(sign(diff( abs(fft_frame).^2) )) == -2) + 1;
        
        % Threshold in Quiet
        A = 3.64 * (f_kHz).^(-0.8) - 6.5 * exp(-0.6 * (f_kHz - 3.3).^2) + (10^(-3))*(f_kHz).^4;
        
        % Masking Spectrum
        k = 1; % factor for downshift
        big_mask = A;
        for peak_count=1:sum(centers * Fs/N <= Fs/2)
            center_freq = centers(peak_count) * Fs/N;
            center_spl = fft_spl(centers(peak_count));
            b = bark(center_freq); % critical band position
            if b >= 14.5 % 14.5 ~ 2.5 KHz
                downshift = k * (14.5 + b);
            else
                downshift = k * (42.5 - b);
            end
            curr_mask = Schroeder(Fs, N, center_freq, center_spl, downshift);
            big_mask = max(big_mask, curr_mask);
        end
        
        %{
        big_mask = A;
        for peak_count=1:sum(noise * Fs/N <= Fs/2)
            center_freq = noise(peak_count) * Fs/N;
            center_spl = fft_spl(noise(peak_count));
            b = bark(center_freq);
            downshift = 14.5 + b;
            curr_mask = Schroeder(Fs, N, center_freq, center_spl, downshift);
            big_mask = max(big_mask, curr_mask);
        end
        
        for peak_count=1:sum(tonal * Fs/N <= Fs/2)
            b = bark(tonal(peak_count) * Fs/N);
            big_mask = max(big_mask,Schroeder(Fs,N,tonal(peak_count)*Fs/N,...
                fft_spl(tonal(peak_count)), 42.5-b));
            center_freq = tonal(peak_count) * Fs/N;
            center_spl = fft_spl(tonal(peak_count));
            b = bark(center_freq);
            downshift = 42.5 - b;
            curr_mask = Schroeder(Fs, N, center_freq, center_spl, downshift);
            big_mask = max(big_mask, curr_mask);
        end
        %}
        
        masks(frame_count, :) = big_mask;
        spls(frame_count, :) = fft_spl;
        
        % allocate later
        % bit_alloc = allocate(New_FFT2,bitrate,scalebits,N,Fs);
        % [Gain,Data] = p_encode(mdct(frames(frame_count,:)),Fs,N,bit_alloc,scalebits);
    end % end of If-Else Statement        
end % end of frame loop

% f_hz = (1:Fs/N:Fs/2);
% tm = 0.008 + 100./f_hz * (0.03 - 0.008); % I don't know how to use this...

for frame_count=length(frames(:, 1)):-1:1
    fft_spl = spls(frame_count, :);
    % simultaneous masking - what we already calculated, how frequencies 
    % mask each other at the same point in time. 
    sim_mask = masks(frame_count, :);
    
    % temporal masking
    temporal_masking = sim_mask;
    for ii=1:N/2
        if frame_count == 1 % no previous frame, not post-masking
            break;
        end
        %c0 = exp(-tm(ii)*((frame_count-1):-1:1));
        c0 = 2.^( -((frame_count-1):-1:1) * N/2/Fs/0.04); % exponential decay
        post_m = max(masks(1:(frame_count - 1), ii)' .* c0); 
        temporal_masking(ii) = max(post_m, temporal_masking(ii));
    end
    
     %big_mask = 20*((temporal_masking/20).^10+(sim_mask/20).^10).^0.1;
     big_mask = max(temporal_masking, sim_mask);
    
    % Signal Spectrum - Masking Spectrum (with max of 0dB)
    New_FFT = fft_spl - big_mask;
    New_FFT_indices = find(New_FFT > 0);
    New_FFT2 = zeros(1, N/2);
    for ii=1:length(New_FFT_indices)
        New_FFT2(New_FFT_indices(ii)) = New_FFT(New_FFT_indices(ii));
    end
    
    if frame_count == -1 % no plots
        semilogx(0:(Fs/2)/(N/2):Fs/2-1,fft_spl(1:N/2),'b');
        hold on;
        semilogx(0:(Fs/2)/(N/2):Fs/2-1,big_mask,'m');
        hold off;
        title('Signal (blue) and Masking Spectrum (pink)');
        figure;
        semilogx(0:(Fs/2)/(N/2):Fs/2-1,New_FFT2);
        title('SMR');
        %figure;
        %stem(allocate(New_FFT2,bitrate,scalebits,N,Fs));
        %title('Bits perceptually allocated');
    end
        
    New_FFT_all(frame_count, :) = New_FFT2;
end
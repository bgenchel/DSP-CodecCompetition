function New_FFT_all = spectrumMasking(frames, N, numBands, Fs, fftmax)
% initialize outoput
New_FFT_all = zeros(length(frames(:,1)), N/2);

for frame_count=1:length(frames(:,1))

    if mod(frame_count, 10) == 0
        outstring = sprintf('Now Masking Frame %i of %i', frame_count, length(frames(:,1)));
        disp(outstring);
    end
    
    fft_frame = fft(frames(frame_count, :));

    if fft_frame == zeros(1, N)
        New_FFT_all(frame_count,:) = zeros(1, N/2);
    else    
        len = length(fft_frame);
        peak_width = zeros(1,len);

        % Find Peaks
        centers = find(diff(sign(diff( abs(fft_frame).^2) )) == -2) + 1;
        spectral_density = zeros(1, length(centers));
    
        peak_max = NaN * ones(size(centers));
        peak_min = NaN*ones(size(centers));
        for k=1:numel(centers) % could also use 'length' here
            peak_max(k) = centers(k) + 2;
            peak_min(k) = centers(k) - 2;
            peak_width(k) = peak_max(k) - peak_min(k); % always 4?
            
            for j=peak_min(k):peak_max(k)
                if (j > 0) && (j < N)
                    spectral_density(k) = spectral_density(k) + abs(fft_frame(j))^2;
                end
            end
        end
        % This gives the squared amplitude of the original signal
        % (this is here just for educational purposes)
        % modified_SD = spectral_density / ((N^2)/8);
        % SPL = 96 + 10 * log10(modified_SD);
        
        % TRANSFORM FFT'S TO SPL VALUES
        fft_spl = 96 + 20 * log10(abs(fft_frame)/fftmax);
    
        % Threshold in Quiet
        f_kHz = (1:Fs/N:Fs/2)/1000;
        A = 3.64 * (f_kHz).^(-0.8) - 6.5 * exp(-0.6 * (f_kHz - 3.3).^2) + (10^(-3))*(f_kHz).^4;
    
        % Masking Spectrum
        big_mask = max(A,Schroeder(Fs,N,centers(1)*Fs/N,fft_spl(centers(1)),...
            14.5 + bark(centers(1)*Fs/N)));
        for peak_count=2:sum(centers * Fs/N<=Fs/2)
            big_mask = max(big_mask,Schroeder(Fs,N,centers(peak_count)*Fs/N,...
                fft_spl(centers(peak_count)), 14.5+bark(centers(peak_count)*Fs/N)));
        end

        % Signal Spectrum - Masking Spectrum (with max of 0dB)
        New_FFT = fft_spl(1:N/2)-big_mask;
        New_FFT_indices = find(New_FFT > 0);
        New_FFT2 = zeros(1,N/2);
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
            figure;
            stem(allocate(New_FFT2,bitrate,scalebits,N,Fs));
            title('Bits perceptually allocated');
        end
        
        New_FFT_all(frame_count, :) = New_FFT2;
        
        % allocate later
        % bit_alloc = allocate(New_FFT2,bitrate,scalebits,N,Fs);
        % [Gain,Data] = p_encode(mdct(frames(frame_count,:)),Fs,N,bit_alloc,scalebits);
    end % end of If-Else Statement        
end % end of frame loop
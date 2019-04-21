function New_mask = Schroeder(Fs, N, freq, spl, downshift)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          SCHROEDER        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the Schroeder masking spectrum for a given frequency and SPL

f_Hz = 1:Fs/N:Fs/2;

% Schroeder Spreading Function
dz = bark(freq) - bark(f_Hz);
mask = 15.81 + 7.5*(dz+0.474) - 17.5*sqrt(1 + (dz+0.474).^2);

New_mask = (mask + spl - downshift);
end
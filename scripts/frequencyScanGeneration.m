function [fs, maxFreq, N1, N2] = frequencyScanGeneration()
% [fs, maxFreq, N1, N2] = frequencyScanGeneration()
% generates a full multisine with a high spectral resolution, to lower the
% SNR. It will be applied at a high sampling frequency to cover up to 1MHz
% and make sure to identify all the poles of the FRF (up to this frequency)
%
% It contains 2 realizations because of how the PXI device works and it
% should have 2 periods to remove the first one (transient removal)
%
% The PXI measurement device seems to have issues with values larger than 1
    desiredRMS = 0.1;
    
    fs = 100e3; % -> nyquist at 50kHz
    maxFreq = 40e3;
    spectralRes1 = 0.25;
    spectralRes2 = 0.1;
    
    % res = fs/N
    N1 = fs/spectralRes1;
    N2 = fs/spectralRes2;
    
        % signal 1
    maxExcBin1 = maxFreq / spectralRes1; % maximum excited frequency bin
    SigFFT = zeros(N1, 1);
    Sig1 = zeros(N1, 2);
    Sel1 = zeros(maxExcBin1, 2);
    
    for rep = 1:2
        SigFFT(2:maxExcBin1+1) = exp(1j * unifrnd(-pi, pi, maxExcBin1, 1));
        Sel1(:, rep) = (1:maxExcBin1)';
    
        Sig1(:, rep) = real(ifft(SigFFT));
    
        % impose the max value
        Sig1(:, rep) = Sig1(:, rep) / rms(Sig1(:, rep)) * desiredRMS;
    end
    
        % signal 2
    maxExcBin2 = maxFreq / spectralRes2; % maximum excited frequency bin
    SigFFT = zeros(N2, 1);
    Sig2 = zeros(N2, 2);
    Sel2 = zeros(maxExcBin2, 2);
    
    for rep = 1:2
        SigFFT(2:maxExcBin2+1) = exp(1j * unifrnd(-pi, pi, maxExcBin2, 1));
        Sel2(:, rep) = (1:maxExcBin2)';
    
        Sig2(:, rep) = real(ifft(SigFFT));
    
        % impose the max value
        Sig2(:, rep) = Sig2(:, rep) / rms(Sig2(:, rep)) * desiredRMS;
    end
    
    
    %% ------- Plot of excitation signal ----------
        figure;
        subplot(2,2,1);
        plot((0:N1-1)/fs, Sig1);
        title('Generated multisine for scanning - time');
        xlabel('Time [s]');
        ylabel('Amplitude [V]');
        subplot(2,2,2);
        plot((0:N1-1)*fs/N1, db(fft(Sig1)), 'o');
        title('Generated multisine for scanning - frequency');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [dBV]');
        subplot(2,2,3);
        plot((0:N2-1)/fs, Sig2);
        title('Generated multisine for scanning - time');
        xlabel('Time [s]');
        ylabel('Amplitude [V]');
        subplot(2,2,4);
        plot((0:N2-1)*fs/N2, db(fft(Sig2)), 'o');
        title('Generated multisine for scanning - frequency');
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [dBV]');
    
    folderPath = '../excitations/';
    signalName = 'frequencyScan';
    save(strcat(folderPath, signalName, '_Sig_E0_S0.mat'), "Sig2");
    save(strcat(folderPath, signalName, '_Sel_E0_S0.mat'), "Sel2");

end
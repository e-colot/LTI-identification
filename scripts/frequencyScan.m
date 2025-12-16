clear; close all; clc;
% generates a full multisine with a high spectral resolution, to lower the
% SNR. It will be applied at a high sampling frequency to cover up to 1MHz
% and make sure to identify all the poles of the FRF (up to this frequency)

% It contains 2 realizations because of how the PXI device works and it
% should have 2 periods to remove the first one (transient removal)

% The PXI measurement device seems to have issues with values larger than 1
desiredMaxVal = 0.99;

fs = 2.2e6; % -> nyquist at 1.1MHz
maxFreq = 1e6;
spectralRes = 0.1;

% res = fs/N
N = fs/spectralRes;

maxExcBin = maxFreq / spectralRes; % maximum excited frequency bin
SigFFT = zeros(N, 1);
Sig = zeros(N, 2);
Sel = zeros(maxExcBin, 2);

for rep = 1:2
    SigFFT(2:maxExcBin+1) = exp(1j * unifrnd(-pi, pi, maxExcBin, 1));
    Sel(:, rep) = (1:maxExcBin)';

    Sig(:, rep) = real(ifft(SigFFT));

    % impose the max value
    Sig(:, rep) = Sig(:, rep) / max(Sig(:, rep)) * desiredMaxVal;
end


%% ------- Plot of excitation signal ----------
    figure;
    subplot(2,1,1);
    plot((0:N-1)/fs, Sig);
    title('Generated multisine for scanning - time');
    xlabel('Time [s]');
    ylabel('Amplitude [V]');
    subplot(2,1,2);
    plot((0:N-1)*fs/N, db(fft(Sig)), 'o');
    title('Generated multisine for scanning - frequency');
    xlabel('Frequency [Hz]');
    ylabel('Amplitude [dBV]');

folderPath = '../excitations/';
signalName = 'frequencyScan';
save(strcat(folderPath, signalName, '_Sig_E0_S0.mat'), "Sig");
save(strcat(folderPath, signalName, '_Sel_E0_S0.mat'), "Sel");

pause;

%% ------- Output display ----------

[u, y, sel, sig, ~, ~] = acquisition("frequencyScan/out");

  % remove the second realization (useless)
    u = u(:, 1);
    y = y(:, 1);
  % transient removal
    assert(length(u) == 2*N);
    u = u(N+1:end);
    y = y(N+1:end);

U = fft(u);    
Y = fft(y);

figure;
subplot(211);
plot((0:N-1)*fs/N, db(U), 'o', LineWidth=2);
xlim([1/fs maxFreq]);
title('Input spectrum');
xlabel('Frequency [Hz]');
ylabel('Amplitude [dBV]');

subplot(212);
plot((0:N-1)*fs/N, db(Y), 'o', LineWidth=2);
xlim([1/fs maxFreq]);
title('Output spectrum');
xlabel('Frequency [Hz]');
ylabel('Amplitude [dBV]');

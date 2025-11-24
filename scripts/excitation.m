clear; close all; clc;

folderPath = '../excitations/';
signalNameOdd = 'oddM';
signalNameFull = 'fullM';

%% RMS + band of interest determination

desiredRMS = 0.1;
N = 25e3;

Nfull = 5e3;

NppOdd = 2; % number of realizations
NppFull = 20;
maxExcBin = N/4; % maximum excited frequency bin

oddSig = zeros(N, NppOdd);
oddSel = [];
fullSig = zeros(Nfull, NppFull);
fullSel = [maxExcBin, NppFull];

% odd
for pp = 1:NppOdd
    if pp == 1
        [sig, sel] = oddMultisine(N, desiredRMS, 'maxExcBin', maxExcBin);
    else
        % keep the same excited frequencies for all realizations
        [sig, ~] = oddMultisine(N, desiredRMS, 'sel', sel, 'maxExcBin', maxExcBin);
    end

    oddSig(:, pp) = sig;
    oddSel(:, pp) = sel;
end

% full
maxExcBin = Nfull/4; % maximum excited frequency bin
fullSigFFT = zeros(Nfull, 1);

for pp = 1:NppFull
    selIndex = 1;
    for bin = 1:maxExcBin
        fullSigFFT(bin + 1) = exp(1j * unifrnd(-pi, pi));
        fullSel(selIndex, pp) = bin;
        selIndex = selIndex + 1;
    end
    fullSig(:, pp) = real(ifft(fullSigFFT));

    % impose the rms
    fullSig(:, pp) = fullSig(:, pp) / rms(fullSig(:, pp)) * desiredRMS;
end


save(strcat(folderPath, signalNameOdd, '_Sig_E0_S0.mat'), "oddSig");
save(strcat(folderPath, signalNameOdd, '_Sel_E0_S0.mat'), "oddSel");

save(strcat(folderPath, signalNameFull, '_Sig_E0_S0.mat'), "fullSig");
save(strcat(folderPath, signalNameFull, '_Sel_E0_S0.mat'), "fullSel");


% DEBUG
figure;
subplot(2,1,1);
plot(oddSig);
title('Generated multisine signal(s) - Odd');
subplot(2,1,2);
stem(abs(fft(oddSig)));
title('FFT of the generated multisine signal(s) - Odd');
figure;
subplot(2,1,1);
plot(fullSig);
title('Generated multisine signal(s) - Full');
subplot(2,1,2);
stem(abs(fft(fullSig)));
title('FFT of the generated multisine signal(s) - Full');

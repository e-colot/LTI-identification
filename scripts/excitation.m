clear; close all; clc;

folderPath = '../excitations/';
signalName = 'multiSine';

%% RMS + band of interest determination

desiredRMS = 0.1;
N = 6e3;

Npp = 40; % number of realizations

totalSig = zeros(N, Npp);
totalSel = [];
maxExcBin = N/4; % maximum excited frequency bin

for pp = 1:Npp
    if pp == 1
        [sig, sel] = oddMultisine(N, desiredRMS, 'maxExcBin', maxExcBin);
    else
        % keep the same excited frequencies for all realizations
        [sig, ~] = oddMultisine(N, desiredRMS, 'sel', sel, 'maxExcBin', maxExcBin);
    end

    totalSig(:, pp) = sig;
    totalSel(:, pp) = sel;
end

save(strcat(folderPath, signalName, '_Sig_E0_S0.mat'), "totalSig");
save(strcat(folderPath, signalName, '_Sel_E0_S0.mat'), "totalSel");

disp(strcat('Signal saved to', folderPath, signalName, '_Sig_E0_S0.mat'));
disp(['With ' num2str(size(totalSig,1)) ' points and ' num2str(Npp) ' realizations.']);
disp(' ');
disp('saved to multiSine_Sel_E0_S0.mat');
disp(strcat('Excited frequencies saved to', folderPath, signalName, '_Sel_E0_S0.mat'));
disp(['With ' num2str(size(totalSel,1)) ' excited frequencies.']);

% DEBUG
figure;
subplot(2,1,1);
plot(totalSig);
title('Generated multisine signal(s)');
subplot(2,1,2);
stem(abs(fft(totalSig)));
title('FFT of the generated multisine signal(s)');

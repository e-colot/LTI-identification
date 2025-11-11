clear; close all; clc;



%% RMS + band of interest determination

desiredRMS = 1;
N = 50e3;

Npp = 2; % number of realizations

totalSig = zeros(N, Npp);
totalSel = [];
maxExcBin = N/20; % maximum excited frequency bin

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

save('multiSine_Sig_E0_S0.mat', "totalSig");
save('multiSine_Sel_E0_S0.mat', "totalSel");

disp('Signal saved to multiSine_Sig_E0_S0.mat');
disp(['With ' num2str(size(totalSig,1)) ' points and ' num2str(Npp) ' realizations.']);
disp(' ');
disp('Excited frequencies saved to multiSine_Sel_E0_S0.mat');
disp(['With ' num2str(size(totalSel,1)) ' excited frequencies.']);

% DEBUG
figure;
subplot(2,1,1);
plot(totalSig);
title('Generated multisine signal(s)');
subplot(2,1,2);
stem(abs(fft(totalSig)));
title('FFT of the generated multisine signal(s)');

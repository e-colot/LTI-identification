clear; close all; clc;

fs = 5e3;

%% Load data
[u, y, realizations, power_levels] = acquisition(fs);

% assuming the excited freqeuencies are the same for all realizations and power levels
load('../excitations/multiSine_Sel_E0_S0.mat'); % loads variable 'totalSel'
load('../excitations/multiSine_Sig_E0_S0.mat'); % loads variable 'totalSig'

N = size(u, 1);
periodN = size(totalSig, 1); % number of samples of the original period
repNumber = N / periodN; % number of periods in the acquired signal

transientPeriods = 1;
assert(transientPeriods < repNumber, 'Transient periods to remove exceed total number of periods.');

u = u(transientPeriods*periodN + 1:end, :, :);
y = y(transientPeriods*periodN + 1:end, :, :);

N = size(u, 1);
repNumber = repNumber - transientPeriods;

excitedFreq = (totalSel(:, 1)*repNumber + 1)';

H_est = zeros(N, realizations, power_levels);
H_est_mean = zeros(N, power_levels);
for p = 1:power_levels
    for r = 1:realizations
        U = fft(u(:, r, p));
        Y = fft(y(:, r, p));
        H_est(:, r, p) = Y ./ U;
    end
    H_est_mean(:, p) = mean(H_est(:, :, p), 2);
end

FRF = H_est_mean(excitedFreq, :);
f = (0:N-1)*fs/N;
f = f(excitedFreq);

figure;
plot(f, db(FRF), 'o', 'LineWidth', 2);
xlim([fs/N max(f)]);
title('Estimated Frequency Response Function (FRF)');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;


%% Parametric model 
% try different models :
% - LS
% - TLS
% - GTLS
% - BTLS
% - ML
% with different orders


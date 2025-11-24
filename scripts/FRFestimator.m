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

%% Remove transient

% transient visualization
figure;
subplot(211);
plot(db(y(1:periodN, 1, 1) - y(periodN + 1:2*periodN, 1, 1)), "LineWidth", 2);
title('y_{P1} - y_{P2}');
xlabel('Samples');
ylabel('Amplitude (dB)');
grid on;

transientPeriods = 1;
assert(transientPeriods < repNumber, 'Transient periods to remove exceed total number of periods.');

u = u(transientPeriods*periodN + 1:end, :, :);
y = y(transientPeriods*periodN + 1:end, :, :);

N = size(u, 1);
repNumber = repNumber - transientPeriods;

%% --------------- FRF estimate ---------------
%% Fast method

[FRF_fast, N_exc, N_even, N_odd, N_noise] = fastMethod(u(:, 1, 1), y(:, 1, 1), totalSel, repNumber);
indices = {N_exc, N_even, N_odd, N_noise};
colors = {'k', 'g', 'r', 'b'};
labels = {'Excited', 'Non-excited Even', 'Non-excited Odd', 'Noise'};

f = (0:N-1)'*(fs/N);

figure; 
for itr = 1:4
    subplot(2, 2, itr);
    plot(f(indices{itr}), db(abs(FRF_fast(indices{itr}))), 'o', 'Color', colors{itr}, 'LineWidth', 2); 
    title(['Fast Method FRF Estimate - ' labels{itr}]);
    xlabel('Frequency (Hz)');
    ylabel('|FRF| (dB)');
    xlim([1/fs fs/2]);
    hold on;
end

%% Robust method

G_2D = zeros(realizations, repNumber, periodN);
G_1D = zeros(realizations, periodN);
noiseVar_2D = zeros(realizations, periodN);

% using a single power level
powerLevel = 1;
% for each realization
for r = 1:realizations
    % for each period
    for p = 1:repNumber
        u_per = u((p-1)*periodN + 1:p*periodN, r, powerLevel);
        y_per = y((p-1)*periodN + 1:p*periodN, r, powerLevel);
        
        U = fft(u_per);
        Y = fft(y_per);

        G_2D(r, p, :) = Y ./ U;
    end
    G_1D(r, :) = squeeze(mean(G_2D(r, :, :), 2));
    noiseVar_2D(r, :) = squeeze(var(G_2D(r, :, :), 0, 2)) / repNumber;
end

% sample mean
G_ML = squeeze(mean(G_1D, 1));
noise_var = mean(noiseVar_2D, 1);
distortion_var = var(G_1D, 0, 1) / realizations;

f = (0:periodN-1)'*(fs/periodN);

figure;
subplot(211);
plot(f, db(G_ML), 'o', 'LineWidth', 2);
title('Robust Method FRF Estimate');
xlabel('Frequency (Hz)');
ylabel('|FRF| (dB)');
xlim([1/fs fs/2]);
grid on;

subplot(212);
hold on;
plot(f, db(noise_var), 'r', 'LineWidth', 2);
plot(f, db(distortion_var), 'b', 'LineWidth', 2);
plot(f, db(noise_var + distortion_var), 'g', 'LineWidth', 2);
title('Variance Estimates');
xlabel('Frequency (Hz)');
ylabel('Variance (dB)');
xlim([1/fs fs/2]);
grid on;
legend('Noise Variance', 'Distortion Variance', 'Total Variance');



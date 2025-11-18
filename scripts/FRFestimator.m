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
transientPeriods = 1;
assert(transientPeriods < repNumber, 'Transient periods to remove exceed total number of periods.');

u = u(transientPeriods*periodN + 1:end, :, :);
y = y(transientPeriods*periodN + 1:end, :, :);

N = size(u, 1);
repNumber = repNumber - transientPeriods;


%% Plot FFTs of input and output
f = (0:N-1)*fs/N;


% Splitting the frequency bins into excited, non-excited even, non-excited odd and noise bins
excitedFreq = (totalSel(:, 1)*repNumber + 1)';
oddFreq = setdiff((1:2:N/repNumber)*repNumber + 1, excitedFreq);
evenFreq = setdiff((0:2:(N-1)/repNumber)*repNumber + 1, excitedFreq);
noiseFreq = setdiff(1:N, [excitedFreq, oddFreq, evenFreq]);

figure;

subplot(2,4,[1 2 5 6]);
U1 = fft(u(:, 1, 1));
plot(f(excitedFreq), db(U1(excitedFreq)), 'o', 'LineWidth', 2);
hold on;
otherFreq = [oddFreq, evenFreq, noiseFreq];
plot(f(otherFreq), db(U1(otherFreq)), 'x', 'LineWidth', 2);
title('Input signal FFT');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('Excited bins', 'Non-excited bins');
xlim([1 fs/2]);
limitYAxis = ylim;
grid on;

subplot(2,4,3);
Y1 = fft(y(:, 1, 1));
plot(f(excitedFreq), db(Y1(excitedFreq)), 'o', 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
title('Output signal FFT, excited frequencies');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
xlim([1 fs/2]);
ylim(limitYAxis);
grid on;
subplot(2,4,4);
plot(f(oddFreq), db(Y1(oddFreq)), 'x', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);
title('Output signal FFT, non-excited odd frequencies');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
xlim([1 fs/2]);
ylim(limitYAxis);
grid on;
subplot(2,4,7);
plot(f(evenFreq), db(Y1(evenFreq)), '*', 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2);
title('Output signal FFT, non-excited even frequencies');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
xlim([1 fs/2]);
ylim(limitYAxis);
grid on;
subplot(2,4,8);
plot(f(noiseFreq), db(Y1(noiseFreq)), '.', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2);
title('Output signal FFT, noise frequencies');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
xlim([1 fs/2]);
ylim(limitYAxis);
grid on;


%% FRF estimate

H_est = zeros(N, realizations, power_levels);
H_est_mean = zeros(N, power_levels);
for p = 1:power_levels
    for r = 1:realizations
        U = fft(u(:, r, p));
        Y = fft(y(:, r, p));
        H_est(:, r, p) = Y ./ U;
    end
    H_est_mean(:, p) = mean(H_est(:, :, p), 2);

    %% Standard deviation on the FRF estimate
    H_var = var(H_est(:, :, p), 0, 2);
    H_std = sqrt(H_var);

    figure("Name", strcat("FRF Estimate, Power Level ", num2str(p-1)));
    % Excited
    subplot(2,2,1);
    plot(f(excitedFreq), db(H_est_mean(excitedFreq, p)), 'o', 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
    hold on;
    fill([f(excitedFreq), fliplr(f(excitedFreq))], ...
         [(db(H_est_mean(excitedFreq, p) + H_std(excitedFreq)))', ...
          fliplr((db(H_est_mean(excitedFreq, p) - H_std(excitedFreq)))')], ...
          'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    title('Excited bins with std dev');
    xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
    xlim([1 fs/2]); ylim(limitYAxis); grid on;

    % Non-excited odd
    subplot(2,2,2);
    plot(f(oddFreq), db(H_est_mean(oddFreq, p)), 'x', 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
    hold on;
    fill([f(oddFreq), fliplr(f(oddFreq))], ...
         [(db(H_est_mean(oddFreq, p) + H_std(oddFreq)))', ...
          fliplr((db(H_est_mean(oddFreq, p) - H_std(oddFreq)))')], ...
          'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    title('Non-excited odd bins with std dev');
    xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
    xlim([1 fs/2]); ylim(limitYAxis); grid on;

    % Non-excited even
    subplot(2,2,3);
    plot(f(evenFreq), db(H_est_mean(evenFreq, p)), '*', 'LineWidth', 2, 'Color', [0.9290 0.6940 0.1250]);
    hold on;
    fill([f(evenFreq), fliplr(f(evenFreq))], ...
         [(db(H_est_mean(evenFreq, p) + H_std(evenFreq)))', ...
          fliplr((db(H_est_mean(evenFreq, p) - H_std(evenFreq)))')], ...
          'y', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    title('Non-excited even bins with std dev');
    xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
    xlim([1 fs/2]); ylim(limitYAxis); grid on;
    
    % Noise
    subplot(2,2,4);
    plot(f(noiseFreq), db(H_est_mean(noiseFreq, p)), '.', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
    hold on;
    fill([f(noiseFreq), fliplr(f(noiseFreq))], ...
         [(db(H_est_mean(noiseFreq, p) + H_std(noiseFreq)))', ...
          fliplr((db(H_est_mean(noiseFreq, p) - H_std(noiseFreq)))')], ...
          [0.5 0.5 0.5], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    title('Noise bins with std dev');
    xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
    xlim([1 fs/2]); ylim(limitYAxis); grid on;

    % Signal PSD
    figure("Name", strcat("Output Signal PSD, Power Level ", num2str(p-1)));
    Syy = (1/(N*fs)) * abs(Y1).^2;
    % Excited
    subplot(2,2,1);
    plot(f(excitedFreq), db(Syy(excitedFreq, p)), 'o', 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
    title('Excited bins');
    xlabel('Frequency (Hz)'); 
    ylabel('Magnitude (dB)');
    xlim([1 fs/2]); 
    ylim([-300, max(db(Syy)) + 10]);
    grid on;

    % Non-excited odd
    subplot(2,2,2);
    plot(f(oddFreq), db(Syy(oddFreq, p)), 'x', 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
    title('Non-excited odd bins');
    xlabel('Frequency (Hz)'); 
    ylabel('Magnitude (dB)');
    xlim([1 fs/2]); 
    ylim([-300, max(db(Syy)) + 10]);
    grid on;

    % Non-excited even
    subplot(2,2,3);
    plot(f(evenFreq), db(Syy(evenFreq, p)), '*', 'LineWidth', 2, 'Color', [0.9290 0.6940 0.1250]);
    title('Non-excited even bins');
    xlabel('Frequency (Hz)'); 
    ylabel('Magnitude (dB)');
    xlim([1 fs/2]); 
    ylim([-300, max(db(Syy)) + 10]);
    grid on;
    
    % Noise
    subplot(2,2,4);
    plot(f(noiseFreq), db(Syy(noiseFreq, p)), '.', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
    title(['Noise bins']);
    xlabel('Frequency (Hz)'); 
    ylabel('Magnitude (dB)');
    xlim([1 fs/2]); 
    ylim([-300, max(db(Syy)) + 10]);
    grid on;

    


end

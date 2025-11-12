clear; close all; clc;

fs = 2e3;
[u, y] = aquisition(fs);

N = size(u, 1);
realizations = size(u, 2);
power_levels = size(u, 3);

f = (0:N-1)*fs/N;

figure;
subplot(2,1,1);
plot(f, db(fft(u(:, 1))));
hold on;
plot(f, db(fft(u(:, 2))));
title('Input signal FFT');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('Realization 1', 'Realization 2');
xlim([0 fs/2]);
grid on;
subplot(2,1,2);
plot(f, db(fft(y(:, 1))));
hold on;
plot(f, db(fft(y(:, 2))));
title('Output signal FFT');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('Realization 1', 'Realization 2');
xlim([0 fs/2]);
grid on;

% estimate FRF
H_est = zeros(N, realizations, power_levels);
H_est_mean = zeros(N, power_levels);
for p = 1:power_levels
    for r = 1:realizations
        U = fft(u(:, r, p));
        Y = fft(y(:, r, p));
        H_est(:, r, p) = Y ./ U;
    end
    H_est_mean(:, p) = mean(H_est(:, :, p), 2);
    
    figure;
    plot(f, db(H_est_mean(:, p)));
    title(['Estimated FRF Magnitude - Power Level ', num2str(p-1)]);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    xlim([0 fs/2]);
    grid on;
end
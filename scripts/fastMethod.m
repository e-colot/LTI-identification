function [] = fastMethod(dataFile, fs)
% fastMethod(dataFile, fs)
%
% Computes the FRF estimate using the fast method.
% Based on a single realization, the fft bins are split in linear contribution, even nonlinear distortion, odd nonlinear distortion and noise bins. 
%

    %% Load data
    [u, y, sel, sig, ~, ~] = acquisition(dataFile);

    N = size(u, 1);
    periodN = size(sig, 1); % number of samples of the original period
    repNumber = N / periodN; % number of periods in the acquired signal

    %% Remove transient
        % visualization
        figure;
        plot(db(y(1:periodN, 1, 1) - y(periodN + 1:2*periodN, 1, 1)), "LineWidth", 2);
        title('y_{P1} - y_{P2}');
        xlabel('Samples');
        ylabel('Amplitude (dB)');
        grid on;

    % transient removal
        transientPeriods = 1;
        assert(transientPeriods < repNumber, 'Transient periods to remove exceed total number of periods.');

        u = u(transientPeriods*periodN + 1:end, :, :);
        y = y(transientPeriods*periodN + 1:end, :, :);

    % update sizes
        N = size(u, 1);
        repNumber = repNumber - transientPeriods;

    %% --------------- FRF estimate ---------------

    % bin categories
        N_exc = (sel(:, 1)*repNumber + 1)';
        N_odd = setdiff((1:2:N/repNumber)*repNumber + 1, N_exc);
        N_even = setdiff((0:2:(N-1)/repNumber)*repNumber + 1, N_exc);
        N_noise = setdiff(1:N, [N_exc, N_odd, N_even]);

    % FRF estimate
        U = fft(u);
        Y = fft(y);

        FRF = Y ./ U;

    %% Plots
        indices = {N_exc, N_even, N_odd, N_noise};
        colors = {'k', 'g', 'r', 'b'};
        labels = {'FRF', 'even distortions', 'odd distortions', 'noise'};

        f = (0:N-1)'*(fs/N);

        figure; 
        for itr = 1:4
            subplot(2, 2, itr);
            plot(f(indices{itr}), db(abs(FRF(indices{itr}))), 'o', 'Color', colors{itr}, 'LineWidth', 2); 
            title(['Fast method - ' labels{itr}]);
            xlabel('Frequency (Hz)');
            ylabel('|FRF| (dB)');
            xlim([1/fs fs/4]);
            hold on;
        end

end
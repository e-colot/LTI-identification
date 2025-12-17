function [] = fastMethod(dataFile, fs, periodsTaken)
% fastMethod(dataFile, fs, [periodsTaken])
%
% Computes the FRF estimate using the fast method.
% Based on a single realization, the fft bins of the output signal spectrum are split in linear contribution, even nonlinear distortion, odd nonlinear distortion and noise bins. 
% periodsTaken allows to limit the number of periods taken

    %% Load data
    [~, y, sel, sig, ~, ~] = acquisition(dataFile);

    N = size(y, 1);
    periodN = size(sig, 1); % number of samples of the original period
    repNumber = N / periodN; % number of periods in the acquired signal

    if (nargin > 2 && repNumber > (periodsTaken+1))
          % periodsTaken+1 because 1 period is removed due to transient
        repNumber = periodsTaken+1;
        N = repNumber * periodN;
        y = y(size(y, 1)-N+1:end, :, :);
    end


    %% Remove transient
        % visualization
        figure;
        plot(db(y(1:periodN, 1, 1) - y(periodN+1:2*periodN, 1, 1)), "LineWidth", 2);
        xlabel('Sample');
        ylabel('Amplitude [dBV]');
        grid on;

    % transient removal
        transientPeriods = 1;
        assert(transientPeriods < repNumber, 'Transient periods to remove exceed total number of periods.');

        y = y(transientPeriods*periodN + 1:end, :, :);

    % update sizes
        N = size(y, 1);
        repNumber = N / periodN; 

    %% --------------- Distinction between types of noise  ---------------

    Y = fft(y);

    % bin categories
        N_exc = (sel(:, 1)*repNumber + 1)';
        N_odd = setdiff((1:2:N/repNumber)*repNumber + 1, N_exc);
        N_even = setdiff((0:2:(N-1)/repNumber)*repNumber + 1, N_exc);
        N_noise = setdiff(1:N, [N_exc, N_odd, N_even]);

    %% Plots
        indices = {N_exc, N_even, N_odd, N_noise};
        colors = {'k', 'g', 'r', 'b'};
        labels = {'Linear', 'Even distortions', 'Odd distortions', 'Noise'};

        f = (0:N-1)'*(fs/N);

        figure; 
        for itr = 1:4
            plot(f(indices{itr}), db(abs(Y(indices{itr}, 1))), 'o', 'Color', colors{itr}, 'LineWidth', 1); 
            hold on;
        end
        
        xlabel('Frequency [Hz]');
        ylabel('|Y| [dBV]');
        xlim([1/fs fs/4]);
        legend(labels);

end
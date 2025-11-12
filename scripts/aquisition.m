function [u, y] = aquisition(fs)
% [u, y] = aquisition(fs)
%
% Load the acquired input and output signals from .mat files for a given
% sampling frequency fs.
% Inputs:
%   fs - Sampling frequency (Hz)
% Outputs:
%   u - Input signal tensor (samples x realizations x power levels)
%   y - Output signal tensor (samples x realizations x power levels)
%
% It assumes the measurements name format is:
% 'out{fs in kHz}k_ACQ_R{realization}_P{power level}_E0_M0_F0.mat'

    resultsPath = '../results/';
    filename = strcat('out', num2str(fs/1000), 'k');
    realizations = 2;
    power_levels = 1;

    % load a single file to get N
    file = strcat(resultsPath, filename, "ACQ_R0_P0_E0_M0_F0.mat");
    load(file);
    N = length(YR0);

    % u and y are tensor with following dimensions:
    % (sample, realization, power level)
    u = zeros(N, realizations, power_levels);
    y = zeros(N, realizations, power_levels);

    % load measured signals
    for r = 0:realizations-1
        for p = 0:power_levels-1
            file = strcat(resultsPath, filename, "ACQ_R", num2str(r), "_P", num2str(p), "_E0_M0_F0.mat");
            load(file);
            u(:, r+1, p+1) = YR0(:);
            y(:, r+1, p+1) = YR1(:);
        end
    end
end




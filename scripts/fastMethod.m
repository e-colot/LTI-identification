function [FRF, N_exc, N_even, N_odd, N_noise] = fastMethod(u, y, totalSel, repNumber)
% [FRF, N_exc, N_even, N_odd, N_noise] = fastMethod(u, y, totalSel, repNumber)
%
% Computes the Frequency Response Function (FRF) estimate using the fast method.
% assumes the transient has been removed.
% Inputs:
%   u          - Input signal (time domain)
%   y          - Output signal (time domain)
%   totalSel   - Matrix containing the selected excited frequency bins
%   repNumber  - Number of repetitions in the multisine signal
%
% Outputs:
%   FRF        - Frequency Response Function estimate
%   N_exc      - Indices of excited frequency bins
%   N_even     - Indices of non-excited even frequency bins
%   N_odd      - Indices of non-excited odd frequency bins
%   N_noise    - Indices of noise frequency bins
%

    % fast method uses a single realization and power level
    u = u(:, 1, 1);
    y = y(:, 1, 1);
    N = length(u);

    N_exc = (totalSel(:, 1)*repNumber + 1)';
    N_odd = setdiff((1:2:N/repNumber)*repNumber + 1, N_exc);
    N_even = setdiff((0:2:(N-1)/repNumber)*repNumber + 1, N_exc);
    N_noise = setdiff(1:N, [N_exc, N_odd, N_even]);

    U = fft(u);
    Y = fft(y);

    FRF = Y ./ U;


end
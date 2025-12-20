clear; close all; clc;

fs = 5e3;

Na = 3;
Nb = 3;

%% Get a non parametric model

[f, G_BLA, total_var] = robustMethod("robustMethod/full_5k", fs);

  % only keep frequencies below 1.25kHz (excited ones)
    valid = f <= 1250;
    valid(1) = 0;
    f = f(valid);
    G_BLA = G_BLA(valid);
    total_var = total_var(valid);

%% Obtain the ML estimate

    [A, B, ~, thetaCov] = ML(G_BLA, total_var, f, Na, Nb);

    poles = roots(A);

    df_dtheta = poles.^(0:Na);
    df_dp = sum(((0:Na).*flipud(A)').*(poles.^(-1:Na-1)), 2);

    dp_dtheta = -df_dtheta./df_dp;

      % split real and imaginary part
      dp_dtheta = [real(dp_dtheta); imag(dp_dtheta)];

    cov_p = dp_dtheta * thetaCov(1:Na+1, 1:Na+1) * dp_dtheta';

      % p_var 1 line for each pole
      % The first column contains the variance of the real part, the second
      % column contains the one the imaginary part.
    p_var = [diag(cov_p(1:Na, 1:Na)), diag(cov_p(Na+1:end, Na+1:end))];

    for p = 1:Na
        disp(['Cov of pole ', num2str(p), ':']);
        disp([num2str(p_var(p, 1)), ' + j', num2str(p_var(p, 2))]);
        disp(' ');
    end


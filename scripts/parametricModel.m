clear; close all; clc;

fs = 5e3;

  % order of the denominator of G
    n_a = 10;
  % order of the numerator of G
    n_b = 10;

%% get a non parametric model

[f, G_ML, total_var] = robustMethod("robustMethod/full_5k", fs);

  % remove NaN's from G_ML (at Nyquist frequency)
    valid = ~isnan(G_ML);
    f = f(valid);
    G_ML = G_ML(valid);
    total_var = total_var(valid);
  % only keep frequencies below 1.25kHz (excited ones)
    valid = f <= 1250;
    valid(1) = 0;
    f = f(valid);
    G_ML = G_ML(valid);
    total_var = total_var(valid);

    Ne = length(G_ML);
    NaList = 1:30;
    NbList = 1:30;

%% Initial estimate using TLS

    TLSest = @(Nb, Na) TLS(G_ML, f, Na, Nb);
    [Na_optTLS, Nb_optTLS, A_optTLS, B_optTLS] = AIC(TLSest, NaList, NbList, Ne);
        % [Atls, Btls, ~] = TLS(G_ML, f, n_a, n_b);

    figure(2);
    subplot(211);
    hold on;
    plotTF(A_optTLS, B_optTLS, f);

%% Initial estimate using GTLS
    
    GTLSest = @(Nb, Na) GTLS(G_ML, total_var, f, Na, Nb);
    [Na_optGTLS, Nb_optGTLS, A_optGTLS, B_optGTLS] = AIC(GTLSest, NaList, NbList, Ne);
        % [Agtls, Bgtls, ~] = GTLS(G_ML, total_var, f, n_a, n_b);

    figure(2);
    subplot(211);
    hold on;
    plotTF(A_optGTLS, B_optGTLS, f);

%% Iterative estimate using BTLS
    theta0 = [flipud(A_optGTLS); flipud(B_optGTLS)];
    itrMax = 100;
    r = 1;
    [Abtls, Bbtls, ~] = BTLS(G_ML, total_var, f, n_a, n_b, itrMax, theta0, r);

    figure(2);
    subplot(211);
    hold on;
    plotTF(Abtls, Bbtls, f);
    legend('Non parametric', ...
        ['TLS estimate, Na=', num2str(Na_optTLS), ' Nb=', num2str(Nb_optTLS)], ...
        ['GTLS estimate, Na=', num2str(Na_optGTLS), ' Nb=', num2str(Nb_optGTLS)], ...
        'BTLS estimate');


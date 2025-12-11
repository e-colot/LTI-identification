clear; close all; clc;

fs = 5e3;

  % order of the denominator of G
    n_a = 5;
  % order of the numerator of G
    n_b = 5;

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


%% Initial estimate using TLS
    [Atls, Btls] = TLS(G_ML, f, n_a, n_b);

    subplot(211);
    hold on;
    plotTF(Atls, Btls, f);

%% Initial estimate using GTLS
    [Agtls, Bgtls] = GTLS(G_ML, total_var, f, n_a, n_b);

    subplot(211);
    hold on;
    plotTF(Agtls, Bgtls, f);

%% Iterative estimate using BTLS
    theta0 = [flipud(Agtls); flipud(Bgtls)];
    itrMax = 20;
    r = 0.5;
    [Abtls, Bbtls] = BTLS(G_ML, total_var, f, n_a, n_b, itrMax, theta0, r);

    figure(2);
    subplot(211);
    hold on;
    plotTF(Abtls, Bbtls, f);
    legend('Non parametric', 'TLS estimate', 'GTLS estimate', 'BTLS estimate');


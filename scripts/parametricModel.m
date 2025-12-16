clear; close all; clc;

fs = 5e3;

  % order of the denominator of G
    n_a = 10;
  % order of the numerator of G
    n_b = 10;

%% get a non parametric model

[f, G_BLA, total_var] = robustMethod("robustMethod/full_5k", fs);

  % only keep frequencies below 1.25kHz (excited ones)
    valid = f <= 1250;
    valid(1) = 0;
    f = f(valid);
    G_BLA = G_BLA(valid);
    total_var = total_var(valid);

    Ne = length(G_BLA);
    NaList = 1:30;
    NbList = 1:30;

%% Initial estimate using LLS

    LLSest = @(Na, Nb) LLS(G_BLA, f, Na, Nb);
    [Na_optLLS, Nb_optLLS, A_optLLS, B_optLLS] = AIC(LLSest, NaList, NbList, Ne);
        % [Alls, Blls, ~] = LLS(G_BLA, f, n_a, n_b);
    
    figure(2);
    subplot(211);
    hold on;
    plotTF(A_optLLS, B_optLLS, f);

%% Initial estimate using TLS

    TLSest = @(Na, Nb) TLS(G_BLA, f, Na, Nb);
    [Na_optTLS, Nb_optTLS, A_optTLS, B_optTLS] = AIC(TLSest, NaList, NbList, Ne);
        % [Atls, Btls, ~] = TLS(G_BLA, f, n_a, n_b);

    figure(2);
    subplot(211);
    hold on;
    plotTF(A_optTLS, B_optTLS, f);

%% Initial estimate using GTLS
    
    GTLSest = @(Na, Nb) GTLS(G_BLA, total_var, f, Na, Nb);
    [Na_optGTLS, Nb_optGTLS, A_optGTLS, B_optGTLS] = AIC(GTLSest, NaList, NbList, Ne);
        % [Agtls, Bgtls, ~] = GTLS(G_BLA, total_var, f, n_a, n_b);

    figure(2);
    subplot(211);
    hold on;
    plotTF(A_optGTLS, B_optGTLS, f);

%% Iterative estimate using BTLS
    itrMax = 100;
    r = 1;
    BTLSest = @(Na, Nb) BTLS(G_BLA, total_var, f, Na, Nb, itrMax, r);
    [Na_optBTLS, Nb_optBTLS, A_optBTLS, B_optBTLS] = AIC(BTLSest, NaList, NbList, Ne);

    figure(2);
    subplot(211);
    hold on;
    plotTF(A_optBTLS, B_optBTLS, f);

    legend('Non parametric', ...
        ['LLS estimate, Na=', num2str(Na_optLLS), ' Nb=', num2str(Nb_optLLS)], ...
        ['TLS estimate, Na=', num2str(Na_optTLS), ' Nb=', num2str(Nb_optTLS)], ...
        ['GTLS estimate, Na=', num2str(Na_optGTLS), ' Nb=', num2str(Nb_optGTLS)], ...
        ['BTLS estimate, Na=', num2str(Na_optBTLS), ' Nb=', num2str(Nb_optBTLS)]);


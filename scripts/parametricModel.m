clear; close all; clc;

fs = 5e3;

%% get a non parametric model

[f, G_BLA, total_var] = robustMethod("robustMethod/full_5k", fs);

  % only keep frequencies below 1.25kHz (excited ones)
    valid = f <= 1250;
    valid(1) = 0;
    f = f(valid);
    G_BLA = G_BLA(valid);
    total_var = total_var(valid);

    Ne = length(G_BLA);
    NaList = 1:15;
    NbList = 1:15;

%% Initial estimate using LLS

    LLSest = @(Na, Nb) LLS(G_BLA, f, Na, Nb);
    [Na_optLLS, Nb_optLLS, A_optLLS, B_optLLS] = AIC(LLSest, NaList, NbList, Ne);
        % [Alls, Blls, ~] = LLS(G_BLA, f, n_a, n_b);

%% Initial estimate using TLS

    TLSest = @(Na, Nb) TLS(G_BLA, f, Na, Nb);
    [Na_optTLS, Nb_optTLS, A_optTLS, B_optTLS] = AIC(TLSest, NaList, NbList, Ne);
        % [Atls, Btls, ~] = TLS(G_BLA, f, n_a, n_b);

%% Initial estimate using GTLS
    
    GTLSest = @(Na, Nb) GTLS(G_BLA, total_var, f, Na, Nb);
    [Na_optGTLS, Nb_optGTLS, A_optGTLS, B_optGTLS] = AIC(GTLSest, NaList, NbList, Ne);
        % [Agtls, Bgtls, ~] = GTLS(G_BLA, total_var, f, n_a, n_b);

%% Iterative estimate using BTLS
    itrMax = 100;
    r = 1;
    BTLSest = @(Na, Nb) BTLS(G_BLA, total_var, f, Na, Nb, itrMax, r);
    [Na_optBTLS, Nb_optBTLS, A_optBTLS, B_optBTLS] = AIC(BTLSest, NaList, NbList, Ne);

%% ML

    MLest = @(Na, Nb) ML(G_BLA, total_var, f, Na, Nb);
    [Na_optML, Nb_optML, A_optML, B_optML] = AIC(MLest, NaList, NbList, Ne);

%% saving data

save("../results/parametricWorkspace.mat");

%% Plots

    % overlay of all the BLA models
    figure;
    subplot(211);
    plot(f, db(G_BLA), 'o', LineWidth=1.5);
    hold on;
    subplot(212);
      % "cheating" to avoid phase jumps between +pi and -pi
    plot(f, angle(G_BLA)+2*pi*(angle(G_BLA) < -1), 'o', LineWidth=1.5);
    hold on;
    plotTF(A_optLLS, B_optLLS, f);
    plotTF(A_optTLS, B_optTLS, f);
    plotTF(A_optGTLS, B_optGTLS, f);
    plotTF(A_optBTLS, B_optBTLS, f);
    plotTF(A_optML, B_optML, f);

    subplot(211);
        legend('Non parametric', 'LLS', 'TLS', 'GTLS', 'BTLS', 'ML');
        xlabel('Frequency [Hz]');
        ylabel('FRF Magnitude [dB]');
        xlim([f(1) f(end)]);
    subplot(212);
        %legend('Non parametric', 'LLS', 'TLS', 'GTLS', 'BTLS', 'ML');
        xlabel('Frequency [Hz]');
        ylabel('FRF Phase [rad]');
        xlim([f(1) f(end)]);


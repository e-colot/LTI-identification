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
    r = 1;
    BTLSest = @(Na, Nb) BTLS(G_BLA, total_var, f, Na, Nb, r);
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

    % model orders on a plot
    figure;
    hold on;
    plot(Na_optLLS, Nb_optLLS, 'o', MarkerSize=8, DisplayName='LLS', LineWidth=2);
    plot(Na_optTLS, Nb_optTLS, 's', MarkerSize=8, DisplayName='TLS', LineWidth=2);
    plot(Na_optGTLS, Nb_optGTLS, '^', MarkerSize=8, DisplayName='GTLS', LineWidth=2);
    plot(Na_optBTLS, Nb_optBTLS, 'd', MarkerSize=8, DisplayName='BTLS', LineWidth=2);
    plot(Na_optML, Nb_optML, 'v', MarkerSize=8, DisplayName='ML', LineWidth=2);
    plot([0 16], [0 16], 'r--','HandleVisibility','off');
    xlabel('Optimal Model Order n_a');
    ylabel('Optimal Model Order n_b');
    legend('Location', 'best');
    grid on;
    xlim([0 16]);
    ylim([0 16]);

        % poles and zeros on plots
    figure;
    hold on;
    poles_LLS = roots(A_optLLS);
    zeros_LLS = roots(B_optLLS);
    colors = lines(5);
    
    plot(real(poles_LLS), imag(poles_LLS), 'x', MarkerSize=8, DisplayName='LLS poles', LineWidth=2, Color=colors(1,:));
    plot(real(zeros_LLS), imag(zeros_LLS), 'o', MarkerSize=8, DisplayName='LLS zeros', LineWidth=2, Color=colors(1,:));

    poles_TLS = roots(A_optTLS);
    zeros_TLS = roots(B_optTLS);
    plot(real(poles_TLS), imag(poles_TLS), 'x', MarkerSize=8, DisplayName='TLS poles', LineWidth=2, Color=colors(2,:));
    plot(real(zeros_TLS), imag(zeros_TLS), 'o', MarkerSize=8, DisplayName='TLS zeros', LineWidth=2, Color=colors(2,:));

    xline(0, 'k--', 'HandleVisibility', 'off');
    xlabel('Real [rad/s]');
    ylabel('Imaginary [rad/s]');
    legend('Location', 'best');
    grid on;
    %title('LLS and TLS');

    figure;
    hold on;
    poles_GTLS = roots(A_optGTLS);
    zeros_GTLS = roots(B_optGTLS);
    plot(real(poles_GTLS), imag(poles_GTLS), 'x', MarkerSize=8, DisplayName='GTLS poles', LineWidth=2, Color=colors(3,:));
    plot(real(zeros_GTLS), imag(zeros_GTLS), 'o', MarkerSize=8, DisplayName='GTLS zeros', LineWidth=2, Color=colors(3,:));

    poles_BTLS = roots(A_optBTLS);
    zeros_BTLS = roots(B_optBTLS);
    plot(real(poles_BTLS), imag(poles_BTLS), 'x', MarkerSize=8, DisplayName='BTLS poles', LineWidth=2, Color=colors(4,:));
    plot(real(zeros_BTLS), imag(zeros_BTLS), 'o', MarkerSize=8, DisplayName='BTLS zeros', LineWidth=2, Color=colors(4,:));

    poles_ML = roots(A_optML);
    zeros_ML = roots(B_optML);
    plot(real(poles_ML), imag(poles_ML), 'x', MarkerSize=8, DisplayName='ML poles', LineWidth=2, Color=colors(5,:));
    plot(real(zeros_ML), imag(zeros_ML), 'o', MarkerSize=8, DisplayName='ML zeros', LineWidth=2, Color=colors(5,:));

    xline(0, 'k--', 'HandleVisibility', 'off');
    xlabel('Real [rad/s]');
    ylabel('Imaginary [rad/s]');
    legend('Location', 'best');
    grid on;
    %title('GTLS, BTLS, and ML');

%% Stability check

    fprintf('Number of poles with real part >= 0:\n');
    fprintf('LLS: %d\n', sum(real(poles_LLS) >= 0));
    fprintf('TLS: %d\n', sum(real(poles_TLS) >= 0));
    fprintf('GTLS: %d\n', sum(real(poles_GTLS) >= 0));
    fprintf('BTLS: %d\n', sum(real(poles_BTLS) >= 0));
    fprintf('ML: %d\n', sum(real(poles_ML) >= 0));

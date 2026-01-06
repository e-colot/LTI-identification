clear; close all; clc;

fs = 5e3;
colors = get(gca, 'colororder');

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
    plot(f, unwrap(angle(G_BLA)), 'o', LineWidth=1.5);
    hold on;
    plotTF(A_optLLS, B_optLLS, f);
    plotTF(A_optTLS, B_optTLS, f);
    plotTF(A_optGTLS, B_optGTLS, f);
    plotTF(A_optBTLS, B_optBTLS, f);
    plotTF(A_optML, B_optML, f);

    subplot(211);
        legend('Non parametric', 'LLS', 'TLS', 'GTLS', 'BTLS', 'ML');
        xlabel('Frequency [Hz]');
        ylabel('|G_{est}| [dB]');
        xlim([f(1) f(end)]);
        grid on;
    subplot(212);
        %legend('Non parametric', 'LLS', 'TLS', 'GTLS', 'BTLS', 'ML');
        xlabel('Frequency [Hz]');
        ylabel('\angle G_{est} [rad]');
        xlim([f(1) f(end)]);
        grid on;

%% model orders on a plot
Na_optLLS = 35; Nb_optLLS = 32; % computed with large Na_max and Nb_max, hardcoded for speed
    figure;
    hold on;
    plot(Na_optLLS, Nb_optLLS, 'o', Color=colors(2, :), MarkerSize=8, DisplayName='LLS', LineWidth=2);
    plot(Na_optTLS, Nb_optTLS, 's', Color=colors(3, :), MarkerSize=8, DisplayName='TLS', LineWidth=2);
    plot(Na_optGTLS, Nb_optGTLS, '^', Color=colors(4, :), MarkerSize=8, DisplayName='GTLS', LineWidth=2);
    plot(Na_optBTLS, Nb_optBTLS, 'd', Color=colors(5, :), MarkerSize=8, DisplayName='BTLS', LineWidth=2);
    plot(Na_optML, Nb_optML, 'v', Color=colors(6, :), MarkerSize=8, DisplayName='ML', LineWidth=2);
    plot([0 40], [0 40], 'r--','HandleVisibility','off');
    xlabel('Optimal Model Order n_a');
    ylabel('Optimal Model Order n_b');
    legend('Location', 'best');
    grid on;
    xlim([0 40]);
    ylim([0 40]);

%% poles and zeros on plots
    figure;
    hold on;
    poles_LLS = roots(A_optLLS);
    zeros_LLS = roots(B_optLLS);
    
    plot(real(poles_LLS), imag(poles_LLS), 'x', Color=colors(2, :), MarkerSize=8, DisplayName='LLS poles', LineWidth=2);
    plot(real(zeros_LLS), imag(zeros_LLS), 'o', Color=colors(2, :), MarkerSize=8, DisplayName='LLS zeros', LineWidth=2);

    poles_TLS = roots(A_optTLS);
    zeros_TLS = roots(B_optTLS);
    plot(real(poles_TLS), imag(poles_TLS), 'x', Color=colors(3, :), MarkerSize=8, DisplayName='TLS poles', LineWidth=2);
    plot(real(zeros_TLS), imag(zeros_TLS), 'o', Color=colors(3, :), MarkerSize=8, DisplayName='TLS zeros', LineWidth=2);

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
    plot(real(poles_GTLS), imag(poles_GTLS), 'x', MarkerSize=8, DisplayName='GTLS poles', LineWidth=2, Color=colors(4, :));
    plot(real(zeros_GTLS), imag(zeros_GTLS), 'o', MarkerSize=8, DisplayName='GTLS zeros', LineWidth=2, Color=colors(4, :));

    poles_BTLS = roots(A_optBTLS);
    zeros_BTLS = roots(B_optBTLS);
    plot(real(poles_BTLS), imag(poles_BTLS), 'x', MarkerSize=8, DisplayName='BTLS poles', LineWidth=2, Color=colors(5,:));
    plot(real(zeros_BTLS), imag(zeros_BTLS), 'o', MarkerSize=8, DisplayName='BTLS zeros', LineWidth=2, Color=colors(5,:));

    poles_ML = roots(A_optML);
    zeros_ML = roots(B_optML);
    plot(real(poles_ML), imag(poles_ML), 'x', MarkerSize=8, DisplayName='ML poles', LineWidth=2, Color=colors(6,:));
    plot(real(zeros_ML), imag(zeros_ML), 'o', MarkerSize=8, DisplayName='ML zeros', LineWidth=2, Color=colors(6,:));

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

%% Compare G_ML with reversed poles

    poles_ML = roots(A_optML);
    stablePoles = -abs(real(poles_ML)) + 1j * imag(poles_ML);
    A_optML_stable = A_optML(1) * poly(stablePoles);
        
    s = 1j*2*pi*f;

        A1_eval = polyval(A_optML, s);
        B_eval = polyval(B_optML, s);
        G1 = B_eval./A1_eval;

        A2_eval = polyval(A_optML_stable, s);
        G2 = B_eval./A2_eval;
        

    figure;
        subplot(211);
        plot(f, db(G1), 'Color', colors(6,:), 'LineWidth', 1.5);
        hold on;
        plot(f, db(G2), 'Color', colors(7,:), 'LineStyle', '--', 'LineWidth', 1.5);
        legend('True ML', 'Stabilized poles');
        xlabel('Frequency [Hz]');
        ylabel('|G_{est}| [dB]');
        xlim([f(1) f(end)]);
        grid on;


        subplot(212);
        plot(f, unwrap(angle(G1)), 'Color', colors(6,:), 'LineWidth', 1.5);
        hold on;
        plot(f, unwrap(angle(G2)), 'Color', colors(7,:), 'LineStyle', '--', 'LineWidth', 1.5);
        %legend('Non parametric', 'LLS', 'TLS', 'GTLS', 'BTLS', 'ML');
        xlabel('Frequency [Hz]');
        ylabel('\angle G_{est} [rad]');
        xlim([f(1) f(end)]);
        grid on;


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

%% Run AIC on ML

    MLest = @(Na, Nb) ML(G_BLA, total_var, f, Na, Nb);
    [Na_optML, Nb_optML, A_optML, B_optML] = AIC(MLest, NaList, NbList, Ne);

%% Run AIC on stabilized ML

    stableMLest = @(Na, Nb) stableRealML(G_BLA, total_var, f, Na, Nb);
    [Na_optStableML, Nb_optStableML, A_optStableML, B_optStableML] = AIC(stableMLest, NaList, NbList, Ne);

%% Plots

    s = 1j*2*pi*f;

        A1_eval = polyval(A_optML, s);
        B1_eval = polyval(B_optML, s);
        G1 = B1_eval./A1_eval;

        A2_eval = prod(s + A_optStableML', 2);
        B2_eval = polyval(B_optStableML, s);
        G2 = B2_eval./A2_eval;
        

    figure;
        subplot(211);
        plot(f, db(G1), 'Color', colors(6,:), 'LineWidth', 1.5);
        hold on;
        plot(f, db(G2), 'Color', colors(7,:), 'LineStyle', '--', 'LineWidth', 1.5);
        legend('True ML', 'Stabilized Poles');
        xlabel('Frequency [Hz]');
        ylabel('FRF Magnitude [dB]');
        xlim([f(1) f(end)]);
        grid on;


        subplot(212);
        plot(f, unwrap(angle(G1)), 'Color', colors(6,:), 'LineWidth', 1.5);
        hold on;
        plot(f, unwrap(angle(G2)), 'Color', colors(7,:), 'LineStyle', '--', 'LineWidth', 1.5);
        %legend('Non parametric', 'LLS', 'TLS', 'GTLS', 'BTLS', 'ML');
        xlabel('Frequency [Hz]');
        ylabel('FRF Phase [rad]');
        xlim([f(1) f(end)]);
        grid on;

%% Position of poles and zeros

    figure;
    hold on;

    rootsB = roots(B_optML);
    polesA = roots(A_optML);
    plot(real(rootsB), imag(rootsB), 'o', LineWidth=2, Color=colors(6,:));
    plot(real(polesA), imag(polesA), 'x', LineWidth=2, Color=colors(6,:));
    
    rootsB = roots(B_optStableML);
    plot(real(rootsB), imag(rootsB), 'o', LineWidth=2, Color=colors(7,:));
    plot(real(A_optStableML), imag(A_optStableML), 'x', LineWidth=2, Color=colors(7,:));
    
    xline(0, 'k--', 'HandleVisibility', 'off');
    xlabel('Real [rad/s]');
    ylabel('Imaginary [rad/s]');
    legend('Zeros ML', 'Poles ML', 'Zeros Stabilized ML', 'Poles Stabilized ML');
    grid on;

    for p = 1:size(A_optStableML, 1)
        disp(['Pole ', num2str(p), ':']);
        disp([num2str(real(A_optStableML(p))), ' + j', num2str(imag(A_optStableML(p)))]);
        disp(' ');
    end



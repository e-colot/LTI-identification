clear; close all; clc;

fs = 5e3;
colors = get(gca, 'colororder');
colors = [colors; 
            [0.0980, 0.0471, 0.8314];
            [0.2157, 0.4667, 0.2078]];

%% get a non parametric model

[f, G_BLA, total_var] = robustMethod("robustMethod/full_5k", fs);

  % only keep frequencies below 1.25kHz (excited ones)
    valid = f <= 1250;
    valid(1) = 0;
    f = f(valid);
    G_BLA = G_BLA(valid);
    total_var = total_var(valid);

    Ne = length(G_BLA);
    NaList = 1:12;
    NbList = 1:12;

%% Run AIC on ML

    MLest = @(Na, Nb) ML(G_BLA, total_var, f, Na, Nb);
    [Na_ML, Nb_ML, A_ML, B_ML] = AIC(MLest, NaList, NbList, Ne);

    [~, ~, costML] = ML(G_BLA, total_var, f, Na_ML, Nb_ML);

%% Run AIC on stable ML

    stableMLest = @(Na, Nb) stableML(G_BLA, total_var, f, Na, Nb);
    [Na_StableML, Nb_StableML, A_StableML, B_StableML] = AIC(stableMLest, NaList, NbList, Ne);

    [~, ~, costStableML] = stableML(G_BLA, total_var, f, Na_StableML, Nb_StableML);

%% Run AIC on stable real ML

    stableRealMLest = @(Na, Nb) stableRealML(G_BLA, total_var, f, Na, Nb);
    [Na_StableRealML, Nb_StableRealML, A_StableRealML, B_StableRealML] = AIC(stableRealMLest, NaList, NbList, Ne);

    [~, ~, costStableRealML] = stableRealML(G_BLA, total_var, f, Na_StableRealML, Nb_StableRealML);

%% Run AIC on optimal ML

    optimalMLest = @(Na, Nb) optimalML(G_BLA, total_var, f, Na, Nb);
    [Na_OptimalML, Nb_OptimalML, A_OptimalML, B_OptimalML] = AIC(optimalMLest, NaList, NbList, Ne);

    [~, ~, costOptimalML] = optimalML(G_BLA, total_var, f, Na_OptimalML, Nb_OptimalML);

%% Plots

    s = 1j*2*pi*f;

    % ML
    A_ML_eval = polyval(A_ML, s);
    B_ML_eval = polyval(B_ML, s);
    G_ML = B_ML_eval./A_ML_eval;

    % Stable ML
    A_StableML_eval = prod(s + A_StableML', 2);
    B_StableML_eval = polyval(B_StableML, s);
    G_StableML = B_StableML_eval./A_StableML_eval;

    % Stable Real ML
    A_StableRealML_eval = prod(s + A_StableRealML', 2);
    B_StableRealML_eval = polyval(B_StableRealML, s);
    G_StableRealML = B_StableRealML_eval./A_StableRealML_eval;

    % Optimal ML
    A_OptimalML_eval = prod(s + A_OptimalML', 2);
    B_OptimalML_eval = polyval(B_OptimalML, s);
    G_OptimalML = B_OptimalML_eval./A_OptimalML_eval;

    % FRF Magnitude and Phase - ML vs Stable ML
    figure;
    subplot(211);
    plot(f, db(G_StableML), 'Color', colors(7,:), 'LineWidth', 1.5);
    hold on;
    plot(f, db(G_ML), '--', 'Color', colors(6,:), 'LineWidth', 1.5);
    legend('Stable ML', 'ML');
    xlabel('Frequency [Hz]');
    ylabel('FRF Magnitude [dB]');
    xlim([f(1) f(end)]);
    grid on;

    subplot(212);
    plot(f, unwrap(angle(G_StableML)), 'Color', colors(7,:), 'LineWidth', 1.5);
    hold on;
    plot(f, unwrap(angle(G_ML)), '--', 'Color', colors(6,:), 'LineWidth', 1.5);
    %legend('Stable ML', 'ML');
    xlabel('Frequency [Hz]');
    ylabel('FRF Phase [rad]');
    xlim([f(1) f(end)]);
    grid on;

    % FRF Magnitude and Phase - ML vs Stable Real ML
    figure;
    subplot(211);
    plot(f, db(G_StableRealML), 'Color', colors(8,:), 'LineWidth', 1.5);
    hold on;
    plot(f, db(G_ML), '--', 'Color', colors(6,:), 'LineWidth', 1.5);
    legend('Stable Real ML', 'ML');
    xlabel('Frequency [Hz]');
    ylabel('FRF Magnitude [dB]');
    xlim([f(1) f(end)]);
    grid on;

    subplot(212);
    plot(f, unwrap(angle(G_StableRealML)), 'Color', colors(8,:), 'LineWidth', 1.5);
    hold on;
    plot(f, unwrap(angle(G_ML)), '--', 'Color', colors(6,:), 'LineWidth', 1.5);
    %legend('Stable Real ML', 'ML');
    xlabel('Frequency [Hz]');
    ylabel('FRF Phase [rad]');
    xlim([f(1) f(end)]);
    grid on;

    % FRF Magnitude and Phase - ML vs Optimal ML
    figure;
    subplot(211);
    plot(f, db(G_OptimalML), 'Color', colors(9,:), 'LineWidth', 1.5);
    hold on;
    plot(f, db(G_ML), '--', 'Color', colors(6,:), 'LineWidth', 1.5);
    legend('Optimal ML', 'ML');
    xlabel('Frequency [Hz]');
    ylabel('FRF Magnitude [dB]');
    xlim([f(1) f(end)]);
    grid on;

    subplot(212);
    plot(f, unwrap(angle(G_OptimalML)), 'Color', colors(9,:), 'LineWidth', 1.5);
    hold on;
    plot(f, unwrap(angle(G_ML)), '--', 'Color', colors(6,:), 'LineWidth', 1.5);
    %legend('Optimal ML', 'ML');
    xlabel('Frequency [Hz]');
    ylabel('FRF Phase [rad]');
    xlim([f(1) f(end)]);
    grid on;

    %% Position of poles and zeros

    figure;
    hold on;

    % ML
    rootsB_ML = roots(B_ML);
    polesA_ML = roots(A_ML);
    plot(real(rootsB_ML), imag(rootsB_ML), 'o', LineWidth=2, Color=colors(6,:));
    plot(real(polesA_ML), imag(polesA_ML), 'x', LineWidth=2, Color=colors(6,:));
    
    % Stable ML
    rootsB_StableML = roots(B_StableML);
    plot(real(rootsB_StableML), imag(rootsB_StableML), 'o', LineWidth=2, Color=colors(7,:));
    plot(real(A_StableML), imag(A_StableML), 'x', LineWidth=2, Color=colors(7,:));
    
    % Stable Real ML
    rootsB_StableRealML = roots(B_StableRealML);
    plot(real(rootsB_StableRealML), imag(rootsB_StableRealML), 'o', LineWidth=2, Color=colors(8,:));
    plot(real(A_StableRealML), imag(A_StableRealML), 'x', LineWidth=2, Color=colors(8,:));
    
    % Optimal ML
    rootsB_OptimalML = roots(B_OptimalML);
    plot(real(rootsB_OptimalML), imag(rootsB_OptimalML), 'o', LineWidth=2, Color=colors(9,:));
    plot(real(A_OptimalML), imag(A_OptimalML), 'x', LineWidth=2, Color=colors(9,:));
    
    xline(0, 'k--', 'HandleVisibility', 'off');
    xlabel('Real [rad/s]');
    ylabel('Imaginary [rad/s]');
    legend('Zeros ML', 'Poles ML', 'Zeros Stable ML', 'Poles Stable ML', 'Zeros Stable Real ML', 'Poles Stable Real ML', 'Zeros Optimal ML', 'Poles Optimal ML');
    grid on;

%% Model orders on a plot
    figure;
    hold on;
    plot(Na_ML, Nb_ML, 'v', Color=colors(6, :), MarkerSize=8, DisplayName='ML', LineWidth=2);
    plot(Na_StableML, Nb_StableML, 's', Color=colors(7, :), MarkerSize=8, DisplayName='Stable ML', LineWidth=2);
    plot(Na_StableRealML, Nb_StableRealML, '^', Color=colors(8, :), MarkerSize=8, DisplayName='Stable Real ML', LineWidth=2);
    plot(Na_OptimalML, Nb_OptimalML, 'd', Color=colors(9, :), MarkerSize=8, DisplayName='Optimal ML', LineWidth=2);
    plot([0 13], [0 13], 'r--','HandleVisibility','off');
    xlabel('Optimal Model Order n_a');
    ylabel('Optimal Model Order n_b');
    legend('Location', 'best');
    grid on;
    xlim([0 13]);
    ylim([0 13]);

    %% costs
    disp(['ML cost: ', num2str(costML)]);
    disp(['Stable ML cost: ', num2str(costStableML)]);
    disp(['Stable Real ML cost: ', num2str(costStableRealML)]);
    disp(['Optimal ML cost: ', num2str(costOptimalML)]);
    
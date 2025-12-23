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
    NaList = 1:8;
    NbList = 1:8;

%% Run AIC on ML

    MLest = @(Na, Nb) ML(G_BLA, total_var, f, Na, Nb);
    [Na_ML, Nb_ML, A_ML, B_ML] = AIC(MLest, NaList, NbList, Ne);

    [~, ~, costML] = ML(G_BLA, total_var, f, Na_ML, Nb_ML);

%% Run AIC on SML

    SMLest = @(Na, Nb) SML(G_BLA, total_var, f, Na, Nb);
    [Na_SML, Nb_SML, A_SML, B_SML] = AIC(SMLest, NaList, NbList, Ne);

    [~, ~, costSML] = SML(G_BLA, total_var, f, Na_SML, Nb_SML);

%% Run AIC on SRML

    SRMLest = @(Na, Nb) SRML(G_BLA, total_var, f, Na, Nb);
    [Na_SRML, Nb_SRML, A_SRML, B_SRML] = AIC(SRMLest, NaList, NbList, Ne);

    [~, ~, costSRML] = SRML(G_BLA, total_var, f, Na_SRML, Nb_SRML);

%% Run AIC on OML

    OMLest = @(Na, Nb) OML(G_BLA, total_var, f, Na, Nb);
    [Na_OML, Nb_OML, A_OML, B_OML] = AIC(OMLest, NaList, NbList, Ne);

    [~, ~, costOML] = OML(G_BLA, total_var, f, Na_OML, Nb_OML);

%% Plots

    s = 1j*2*pi*f;

    % ML
    A_ML_eval = polyval(A_ML, s);
    B_ML_eval = polyval(B_ML, s);
    G_ML = B_ML_eval./A_ML_eval;

    % Stable ML
    A_SML_eval = prod(s + A_SML', 2);
    B_SML_eval = polyval(B_SML, s);
    G_SML = B_SML_eval./A_SML_eval;

    % Stable Real ML
    A_SRML_eval = prod(s + A_SRML', 2);
    B_SRML_eval = polyval(B_SRML, s);
    G_SRML = B_SRML_eval./A_SRML_eval;

    % Optimal ML
    A_OML_eval = prod(s + A_OML', 2);
    B_OML_eval = polyval(B_OML, s);
    G_OML = B_OML_eval./A_OML_eval;

    % FRF Magnitude and Phase - ML vs SML
    figure;
    subplot(211);
    plot(f, db(G_SML), 'Color', colors(7,:), 'LineWidth', 1.5);
    hold on;
    plot(f, db(G_ML), '--', 'Color', colors(6,:), 'LineWidth', 1.5);
    legend('SML', 'ML');
    xlabel('Frequency [Hz]');
    ylabel('FRF Magnitude [dB]');
    xlim([f(1) f(end)]);
    grid on;

    subplot(212);
    plot(f, unwrap(angle(G_SML)), 'Color', colors(7,:), 'LineWidth', 1.5);
    hold on;
    plot(f, unwrap(angle(G_ML)), '--', 'Color', colors(6,:), 'LineWidth', 1.5);
    %legend('Stable ML', 'ML');
    xlabel('Frequency [Hz]');
    ylabel('FRF Phase [rad]');
    xlim([f(1) f(end)]);
    grid on;

    % FRF Magnitude and Phase - ML vs SRML
    figure;
    subplot(211);
    plot(f, db(G_SRML), 'Color', colors(8,:), 'LineWidth', 1.5);
    hold on;
    plot(f, db(G_ML), '--', 'Color', colors(6,:), 'LineWidth', 1.5);
    legend('SRML', 'ML');
    xlabel('Frequency [Hz]');
    ylabel('FRF Magnitude [dB]');
    xlim([f(1) f(end)]);
    grid on;

    subplot(212);
    plot(f, unwrap(angle(G_SRML)), 'Color', colors(8,:), 'LineWidth', 1.5);
    hold on;
    plot(f, unwrap(angle(G_ML)), '--', 'Color', colors(6,:), 'LineWidth', 1.5);
    %legend('Stable Real ML', 'ML');
    xlabel('Frequency [Hz]');
    ylabel('FRF Phase [rad]');
    xlim([f(1) f(end)]);
    grid on;

    % FRF Magnitude and Phase - ML vs OML
    figure;
    subplot(211);
    plot(f, db(G_OML), 'Color', colors(9,:), 'LineWidth', 1.5);
    hold on;
    plot(f, db(G_ML), '--', 'Color', colors(6,:), 'LineWidth', 1.5);
    legend('OML', 'ML');
    xlabel('Frequency [Hz]');
    ylabel('FRF Magnitude [dB]');
    xlim([f(1) f(end)]);
    grid on;

    subplot(212);
    plot(f, unwrap(angle(G_OML)), 'Color', colors(9,:), 'LineWidth', 1.5);
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
    rootsB_SML = roots(B_SML);
    plot(real(rootsB_SML), imag(rootsB_SML), 'o', LineWidth=2, Color=colors(7,:));
    plot(real(A_SML), imag(A_SML), 'x', LineWidth=2, Color=colors(7,:));
    
    % Stable Real ML
    rootsB_SRML = roots(B_SRML);
    plot(real(rootsB_SRML), imag(rootsB_SRML), 'o', LineWidth=2, Color=colors(8,:));
    plot(real(A_SRML), imag(A_SRML), 'x', LineWidth=2, Color=colors(8,:));
    
    % Optimal ML
    rootsB_OML = roots(B_OML);
    plot(real(rootsB_OML), imag(rootsB_OML), 'o', LineWidth=2, Color=colors(9,:));
    plot(real(A_OML), imag(A_OML), 'x', LineWidth=2, Color=colors(9,:));
    
    xline(0, 'k--', 'HandleVisibility', 'off');
    xlabel('Real [rad/s]');
    ylabel('Imaginary [rad/s]');
    legend('Zeros ML', 'Poles ML', 'Zeros SML', 'Poles SML', 'Zeros SRML', 'Poles SRML', 'Zeros OML', 'Poles OML');
    grid on;

%% Model orders on a plot
    figure;
    hold on;
    plot(Na_ML, Nb_ML, 'v', Color=colors(6, :), MarkerSize=8, DisplayName='ML', LineWidth=2);
    plot(Na_SML, Nb_SML, 's', Color=colors(7, :), MarkerSize=8, DisplayName='SML', LineWidth=2);
    plot(Na_SRML, Nb_SRML, '^', Color=colors(8, :), MarkerSize=8, DisplayName='SRML', LineWidth=2);
    plot(Na_OML, Nb_OML, 'd', Color=colors(9, :), MarkerSize=8, DisplayName='OML', LineWidth=2);
    plot([0 13], [0 13], 'r--','HandleVisibility','off');
    xlabel('Optimal Model Order n_a');
    ylabel('Optimal Model Order n_b');
    legend('Location', 'best');
    grid on;
    xlim([0 6]);
    ylim([0 6]);

    %% costs
    disp(['ML cost: ', num2str(costML)]);
    disp(['SML cost: ', num2str(costSML)]);
    disp(['SRML cost: ', num2str(costSRML)]);
    disp(['OML cost: ', num2str(costOML)]);
    
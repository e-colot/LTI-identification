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

%% Run AIC on stabilized ML

    MLest = @(Na, Nb) stableML(G_BLA, total_var, f, Na, Nb);
    [Na_optML, Nb_optML, A_optML, B_optML] = AIC(MLest, NaList, NbList, Ne);

%% Plots

    figure;
    subplot(211);
    hold on;
    subplot(212);
    hold on;
    plotTF(A_optML, B_optML, f, 'mult');
    [A, B, ~] = ML(G_BLA, total_var, f, 3, 3);
    plotTF(A, B, f, 'poly');

    subplot(211);
        legend('ML', 'ML with stable poles');
        xlabel('Frequency [Hz]');
        ylabel('FRF Magnitude [dB]');
        xlim([f(1) f(end)]);
        grid on;
    subplot(212);
        xlabel('Frequency [Hz]');
        ylabel('FRF Phase [rad]');
        xlim([f(1) f(end)]);
        grid on;



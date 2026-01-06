clear; %close all; clc;

fs = 5e3;

%% Get a non parametric model

[f, G_BLA, total_var] = robustMethod("robustMethod/full_5k", fs);

  % only keep frequencies below 1.25kHz (excited ones)
    valid = f <= 1250;
    valid(1) = 0;
    f = f(valid);
    G_BLA = G_BLA(valid);
    total_var = total_var(valid);

%% ML uncertainty estimate
    [A_ML, ~, ~, p_var_ML] = ML(G_BLA, total_var, f, 3, 3);

%% stable real ML uncertainty estimate
    [A_SRML, ~, ~, p_var_SRML] = SRML(G_BLA, total_var, f, 6, 5);

%% optimal ML uncertainty estimate
    [A_OML, ~, ~, p_var_OML] = OML(G_BLA, total_var, f, 6, 5);

%% Disps

    poles = roots(A_ML);
    disp('---------- ML estimator ----------');
    for p = 1:length(poles)
        disp(['Standard deviation of pole ', num2str(p), ' in ', num2str(poles(p)), ':']);
        disp(num2str(sqrt(p_var_ML(p, 1)) + 1j *sqrt(p_var_ML(p, 2))));
        disp(' ');
    end
    disp(' ');disp(' ');


    poles = A_SRML;
    disp('---------- SRML estimator ----------');
    for p = 1:length(poles)
        disp(['Standard deviation of pole ', num2str(p), ' in ', num2str(poles(p)), ':']);
        disp(num2str(sqrt(p_var_SRML(p, 1)) + 1j *sqrt(p_var_SRML(p, 2))));
        disp(' ');
    end
    disp(' ');disp(' ');


    poles = A_OML;
    disp('---------- OML estimator ----------');
    for p = 1:length(poles)
        disp(['Standard deviation of pole ', num2str(p), ' in ', num2str(poles(p)), ':']);
        disp(num2str(sqrt(p_var_OML(p, 1)) + 1j *sqrt(p_var_OML(p, 2))));
        disp(' ');
    end
    disp(' ');disp(' ');

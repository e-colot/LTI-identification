clear; close all; clc;

resultsPath = '../results/';
filename = 'out2k';
realizations = 1;
power_levels = 1;

% load reference signals
for r = 0:realizations-1
    for p = 0:power_levels-1
        file = strcat(resultsPath, "REF", filename, "ACQ_R", num2str(r), "_P", num2str(p), "_E0_M0_F0.mat");
        load(file);
    end
end


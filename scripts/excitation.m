clear; close all; clc;



%% RMS + band of interest determination

N = 50e3;

randInput = randn(N,1);
randInput = randInput / rms(randInput);   % impose an rms of 1V

save('../excitation signals/randInput.mat','randInput');





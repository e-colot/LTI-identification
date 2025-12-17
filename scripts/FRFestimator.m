clear; close all; clc;

fs = 5e3;

%% fast method

fastMethod("fastMethod1/odd_5k", fs);
fastMethod("fastMethodEven/even", fs, 2);

%% Robust method 

robustMethod("robustMethod/full_5k", fs);


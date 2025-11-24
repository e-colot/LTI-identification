%clear; close all; clc;

fs = 5e3;

%% fast method

fastMethod("fastMethod1/odd_5k", fs);

%% Robust method

robustMethod("robustMethod/full_5k", fs);


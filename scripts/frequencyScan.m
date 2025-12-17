clear; close all; clc;

[fs, maxFreq, N1, N2] = frequencyScanGeneration();

[u1, y1, ~, ~, ~, ~] = acquisition("frequencyScan/test");
[u2, y2, ~, ~, ~, ~] = acquisition("frequencyScan/test2");

uCat = {u1, u2};
yCat = {y1, y2};
NCat = {N1, N2};

for sig = 1:2
    u = uCat{sig};
    y = yCat{sig};
    N = NCat{sig};
      % remove the second realization (useless)
        u = u(:, 1);
        y = y(:, 1);
      % transient removal
        assert(length(u) == 2*N);
        u = u(N+1:end);
        y = y(N+1:end);
    
    U = fft(u);    
    Y = fft(y);
    
    figure;
    subplot(211);
    plot((0:N-1)*fs/(1000*N), db(U), 'o', Color=[0.4940, 0.1840, 0.5560], LineWidth=2);
    xlim(1e-3*[1/fs maxFreq]);
    %title('Input spectrum');
    xlabel('Frequency [kHz]');
    ylabel('Amplitude [dBV]');
    
    subplot(212);
    plot((0:N-1)*fs/(1000*N), db(Y), 'o', Color=[0.4940, 0.1840, 0.5560], LineWidth=2);
    xlim(1e-3*[1/fs maxFreq]);
    %title('Output spectrum');
    xlabel('Frequency [kHz]');
    ylabel('Amplitude [dBV]');
end


%% Time domain plots

    figure;
    subplot(121);
    plot((0:N2-1)/fs, u2(N2+1:end, 1), Color=[0.4660, 0.6740, 0.1880], LineWidth=2);
    xlabel('Time [s]');
    ylabel('Amplitude [V]');
    subplot(122);
    plot((0:N2-1)/fs, y2(N2+1:end, 1), Color=[0.4660, 0.6740, 0.1880], LineWidth=2);
    xlabel('Time [s]');
    ylabel('Amplitude [V]');

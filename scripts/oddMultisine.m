function [sig, sel] = oddMultisine(N, desiredRMS, varargin)
% [sig, sel] = oddMultisine(N, desiredRMS, [maxExcBin])
% Generates an odd multisine signal of length N with desired RMS value.
%
% sig is the time-domain signal
% sel is the vector of excited frequency bins
%
% If maxExcBin is provided, the multisine will only excite frequencies up
% to that bin (default is N/4).
% If sel is provided, only the frequencies in sel will be excited.
%
% 1 out of every 4 odd bins is randomly not excited to measure odd nonlinear effects.
    
    % parse optional arguments: allow (maxExcBin) or (sel) or (maxExcBin, sel)
    selProvided = false;
    maxExcBin = N/4;

    for i = 1:2:length(varargin)
        if strcmp(varargin{i}, 'maxExcBin')
            maxExcBin = varargin{i+1};
        elseif strcmp(varargin{i}, 'sel')
            sel = varargin{i+1}(:);
            selProvided = true;
        end
    end

    % if a selection of bins was provided, generate only those and return
    if selProvided
        multisine_f = zeros(N, 1);

        % sanitize sel: integers, within range, keep odd bins only
        sel = round(sel(:));
        sel = sel(sel >= 1 & sel <= maxExcBin & mod(sel,2) == 1);

        for k = 1:numel(sel)
            bin = sel(k);
            multisine_f(bin + 1) = exp(1j * unifrnd(-pi, pi));
        end

        sig = real(ifft(multisine_f));
        sig = sig * desiredRMS / rms(sig);
        return;
    end

    multisine_f = zeros(N, 1);

    % selected frequencies
    sel = NaN(maxExcBin, 1);
    selIndex = 1;

    % excite only the odd frequencies and remove 1 bin in every group of 4
    for i = 1:maxExcBin/(2*4)
        % choose which bin to remove in the batch
        toSkip = randi([0 3]);
        for j = 0:3
            bin = 2*(4*(i-1) + j) + 1;
            if (j ~= toSkip)
                multisine_f(bin + 1) = exp(1j * unifrnd(-pi, pi));
                sel(selIndex) = bin;
                selIndex = selIndex + 1;
            end
        end
    end

    sel = sel(~isnan(sel));

    sig = real(ifft(multisine_f));

    sig = sig * desiredRMS / rms(sig);

end
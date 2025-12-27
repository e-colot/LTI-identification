function [A, B, cost, costHandle] = LLS(G_BLA, f, Na, Nb)
% [A, B, cost, costHandle] = LLS(G_BLA, f, Na, Nb)
% Determines the parameters A and B using a LLS estimator (with a_0 = 1).
% 
% Na and Nb are the order of the denominator and numerator of the
% parametric transfer function
%
% Written by E. Colot on Dec 15 2025

    % solves y = H * theta
    % H = [sY ... U sU ...]
    % taking U = 1 and Y = G_BLA
    
      % placing s in the matrix H
        H0 = repmat(1j*2*pi*f, 1, Na+Nb+1);
        H0 = H0.^([(1:Na) (0:Nb)]);
      % placing G_BLA in the matrix
        H0 = H0.*[repmat(-G_BLA, 1, Na), -1*ones(size(G_BLA, 1), Nb+1)];
    
      % To force real parameters, the real and imaginary part of J must be
      % split
        H0 = [real(H0)  ;
             imag(H0) ];
    
      % for conditioning, rms normalization of H (column-wise)
        S = diag(1./rms(H0));
        H = H0 * S;

        y = [real(G_BLA); imag(G_BLA)];

        theta_LLS = H\y;

        % because of the rms normalization
        theta_LLS = S * theta_LLS;
        
        A = flipud([1; theta_LLS(1:Na)]);
        B = flipud(theta_LLS(Na+1:end));

%% Cost computation

    cost = sum(abs(H0 * theta_LLS - y).^2);

    costHandle = @(G_BLA, useless, f) costComputation(G_BLA, useless, f);

    function cost = costComputation(G_BLA, ~, f)
        H0 = repmat(1j*2*pi*f, 1, Na+Nb+1);
        H0 = H0.^([(1:Na) (0:Nb)]);
        H0 = H0.*[repmat(-G_BLA, 1, Na), -1*ones(size(G_BLA, 1), Nb+1)];
        H0 = [real(H0)  ;
             imag(H0) ];
        y = [real(G_BLA); imag(G_BLA)];
        cost = sum(abs(H0 * theta_LLS - y).^2);
    end

end

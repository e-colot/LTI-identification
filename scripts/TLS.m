function [A, B, cost, costHandle] = TLS(G_BLA, f, Na, Nb)
% [A, B, cost, costHandle] = TLS(G_BLA, f, Na, Nb)
% Determines the parameters A and B using a TLS estimator.
% 
% Na and Nb are the order of the denominator and numerator of the
% parametric transfer function
%
% Written by E. Colot on Dec 10 2025


    % J = [Y sY ... U sU ...]
    % taking U = 1 and Y = G_BLA
    
      % placing s in the matrix J
        J = repmat(1j*2*pi*f, 1, Na+Nb+2);
        J = J.^([(0:Na) (0:Nb)]);
      % placing G_BLA in the matrix
        J = J.*[repmat(G_BLA, 1, Na+1), -1*ones(size(G_BLA, 1), Nb+1)];
    
      % To force real parameters, the real and imaginary part of J must be
      % split
        J = [real(J)  ;
             imag(J) ];
    
      % for conditioning, rms normalization of J (column-wise)
        S = diag(1./rms(J));
        J = J * S;

        [~,~,V] = svd(J);
        theta_TLS = V(:, end);

        % because of the rms normalization
        theta_TLS = S * theta_TLS;
        
        A = flipud(theta_TLS(1:Na+1));
        B = flipud(theta_TLS(Na+2:end));

%% Cost computation

    s = 1j*2*pi*f;
    A_eval = polyval(A, s);
    B_eval = polyval(B, s);
    
    G_est = B_eval./A_eval;
    err = G_BLA-G_est;

    cost = sum(abs(err).^2);

    costHandle = @(G_BLA, useless, f) costComputation(G_BLA, useless, f);

    function cost = costComputation(G_BLA, ~, f)
        s = 1j*2*pi*f;
        A_eval = polyval(A, s);
        B_eval = polyval(B, s);
        
        G_est = B_eval./A_eval;
        err = G_BLA-G_est;
    
        cost = sum(abs(err).^2);     
    end

end

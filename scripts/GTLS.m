function [A, B, cost] = GTLS(G_BLA, G_var, f, Na, Nb)
% [A, B, cost] = GTLS(G_BLA, G_var, f, Na, Nb)
% Determines the parameters A and B using a GTLS estimator.
% 
% Na and Nb are the order of the denominator and numerator of the
% parametric transfer function
%
% Written by E. Colot on Dec 11 2025

%% ------------ Building J ---------------------

    % J = [Y sY ... U sU ...]
    % taking U = 1 and Y = G_BLA
    
      % placing s in the matrix J
        J = repmat(1j*2*pi*f, 1, Na+Nb+2);
        J = J.^([(0:Na) (0:Nb)]);
      % placing G_BLA in the matrix
        J = J.*[repmat(G_BLA, 1, Na+1), -ones(size(G_BLA, 1), Nb+1)];
    
      % To force real parameters, the real and imaginary part of J must be
      % split
        J = [real(J)  ;
             imag(J) ];
    
      % for conditioning, rms normalization of J (column-wise)
        S = diag(1./rms(J));
        J = J * S;

%% ------------ Building C_J^0.5 ---------------------

    % C_J = column covariance of DeltaJ
        deltaH = zeros(length(f), Na+Nb+2);
        deltaH(:, 1:Na+1) = repmat(-1j*2*pi*f, 1, Na+1).^(0:Na);
        deltaH = deltaH * S;
        deltaH = deltaH .* repmat(sqrt(G_var), 1, Na+Nb+2);
    
        C_J = deltaH' * deltaH;
        
      % To impose real parameters
        C_J = real(C_J);
    
        % To compute C_J^0.5, it must be diagonalized
        [eigenVect, eigenVal] = eig(C_J);
        sqrtC_J = eigenVect * sqrt(eigenVal) * eigenVect';

%% ------------- Computation of theta -------------

        [~,~,XJ,CJ,SJ] = gsvd(J, sqrtC_J);
        
        [~, i] = min(diag(CJ)./diag(SJ));
        
        Xinv = inv(XJ');
        theta_GTLS = 1/SJ(i, i) * Xinv(:, i);
        
        % because of the rms normalization
        theta_GTLS = S * theta_GTLS;
        
        A = flipud(theta_GTLS(1:Na+1));
        B = flipud(theta_GTLS(Na+2:end));

%% Cost computation

    err = J*theta_GTLS;
    cost = sum(abs(err).^2);

end

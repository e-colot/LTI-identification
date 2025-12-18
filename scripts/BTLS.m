function [A, B, cost] = BTLS(G_BLA, G_var, f, Na, Nb, itrMax, r)
% [A, B, cost] = BTLS(G_BLA, G_var, f, Na, Nb, itrMax, r)
% Determines the parameters A and B using a BTLS estimator. The starting
% value for theta is computed using a GTLS estimator
% 
% Na and Nb are the order of the denominator and numerator of the
% parametric transfer function
%
% Written by E. Colot on Dec 13 2025

[Agtls, Bgtls, ~] = GTLS(G_BLA, G_var, f, Na, Nb);
theta0 = [flipud(Agtls); flipud(Bgtls)];

prevCost = NaN;
s = 1j*2*pi*f;

% J0 = [Y sY ... U sU ...]
% taking U = 1 and Y = G_BLA

  % placing s in the matrix J
    J0 = repmat(s, 1, Na+Nb+2);
    J0 = J0.^([(0:Na) (0:Nb)]);
  % placing G_BLA in the matrix
    J0 = J0.*[repmat(G_BLA, 1, Na+1), -ones(size(G_BLA, 1), Nb+1)];


% C_J = column covariance of DeltaJ
    DeltaJ0 = zeros(length(f), Na+Nb+2);
    DeltaJ0(:, 1:Na+1) = repmat(s, 1, Na+1).^(0:Na);

    for itr = 1:itrMax

%% -------- Computation of the weight W ---------
        oldTheta = theta0;

        Aprev = polyval(flipud(oldTheta(1:Na+1)), s);
        Aprev = diag(Aprev);

        var_e = abs(Aprev).^2 * G_var;
        W = diag(var_e.^(-r));

%% ------------ Building J ---------------------

      % Add weight to J
        J = W * J0;
    
      % To force real parameters, the real and imaginary part of J must be
      % split
        J = [real(J)  ;
             imag(J) ];
    
      % for conditioning, rms normalization of J (column-wise)
        S = diag(1./rms(J));
        J = J * S;

%% ------------ Building C_J^0.5 ---------------------

        DeltaJ = DeltaJ0 * S;
        DeltaJ = W * DeltaJ;
        DeltaJ = DeltaJ .* repmat(sqrt(G_var), 1, Na+Nb+2);
    
        C_J = DeltaJ' * DeltaJ;
        
      % To impose real parameters
        C_J = real(C_J);
    
        % To compute C_J^0.5, it must be diagonalized
        [eigenVect, eigenVal] = eig(C_J);
        sqrtC_J = eigenVect * sqrt(eigenVal) * eigenVect';

%% ------------- Computation of theta -------------

        [~,~,XJ,CJ,SJ] = gsvd(J, sqrtC_J);
        
        [~, i] = min(diag(CJ)./diag(SJ));
        
        Xinv = inv(XJ');
        newTheta = 1/SJ(i, i) * Xinv(:, i);
        
        % because of the rms normalization
        newTheta = S * newTheta;
        
%% Check the new cost
        A = flipud(newTheta(1:Na+1));
        B = flipud(newTheta(Na+2:end));
        A_eval = polyval(A, s);
        B_eval = polyval(B, s);
        
        G_est = B_eval./A_eval;
        err = G_BLA-G_est;

        cost = sum(abs(err).^2 ./ (var_e.^r));
        if (cost >= prevCost)
            disp(['BTLS stopped after ', num2str(itr), ' iterations due to convergence']);
            return
        end
        prevCost = cost;

    end

end

function [A, B, cost] = BTLS(G_BLA, G_var, f, Na, Nb, r)
% [A, B, cost] = BTLS(G_BLA, G_var, f, Na, Nb, r)
% Determines the parameters A and B using a BTLS estimator. The starting
% value for theta is computed using a GTLS estimator
% 
% Na and Nb are the order of the denominator and numerator of the
% parametric transfer function
%
% Written by E. Colot on Dec 13 2025

[Agtls, Bgtls, ~] = GTLS(G_BLA, G_var, f, Na, Nb);
oldTheta = [flipud(Agtls); flipud(Bgtls)];

oldCost = Inf;
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
    DeltaJ0 = DeltaJ0 .* repmat(sqrt(G_var), 1, Na+Nb+2);

    itrCnt = 1;

    while 1
    
        %% ------ Weights ------
            % var_e = |A|^2 * G_var
            A = s.^(0:Na)*oldTheta(1:Na+1);
            var_e = abs(A).^2 .* G_var;

            W = diag(var_e.^(-r));
         
        %% ------ Scaling ------
            % for J
            J = W * J0;
            J_re = [real(J); imag(J)];

            S = diag(1./rms(J_re));
            J_re = J_re * S;

            % for C_J
            DeltaJ = W * DeltaJ0 * S;

            C_J_re = real(DeltaJ' * DeltaJ);
                    
            % To compute C_J^0.5, it must be diagonalized
            [eigenVect, eigenVal] = eig(C_J_re);
            sqrt_C_J_re = eigenVect * sqrt(eigenVal) * eigenVect';

        %% ------ Computation of theta ------
            [~,~,XJ,CJ,SJ] = gsvd(J_re, sqrt_C_J_re);
            
            [~, i] = min(diag(CJ)./diag(SJ));
            
            Xinv = inv(XJ');
            newTheta = 1/SJ(i, i) * Xinv(:, i);
            
            % because of the rms normalization
            newTheta = S * newTheta;

        %% ------ Cost computation ------
            A = s.^(0:Na)*newTheta(1:Na+1);
            B = s.^(0:Nb)*newTheta(Na+2:end);

            e = A.*G_BLA - B;
            var_e = abs(A).^2 .* G_var;

            newCost = sum(abs(e).^2 ./ (var_e.^r));

        if (newCost < oldCost)
            oldTheta = newTheta;
            oldCost = newCost;
            itrCnt = itrCnt + 1;
        else
            % disp(['BTLS stopped after ', num2str(itrCnt), ' iterations.']);
            break;
        end
    end

    A = flipud(oldTheta(1:Na+1));
    B = flipud(oldTheta(Na+2:end));
    cost = oldCost;

end
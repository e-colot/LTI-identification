function [A, B, cost, p_var] = OML(G_BLA, G_var, f, Na, Nb, show)
% [A, B, cost, p_var] = SML(G_BLA, G_var, f, Na, Nb, [show])
% Determines the parameters A and B using a ML estimator that imposes 
% stable poles. The starting value for theta is computed using a GTLS 
% estimator
% 
% Na and Nb are the order of the denominator and numerator of the
% parametric transfer function
%
% show is 0 by default and shows how theta changes with the iterations
% p_var returns the covariance of roots(A)
%
% Written by E. Colot on Dec 21 2025

    if nargin <= 5
        show = 0;
    end

    s = 1j*2*pi*f;
    [oldA, oldB, ~] = SRML(G_BLA, G_var, f, Na, Nb);
    oldB = flipud(oldB);

    itrCnt = 1;

    thetaChange = 0;

    % split real and complex poles
    complexPoles = oldA(imag(oldA) > 0);
      % complex ones
        oldAlpha = sqrt(-real(complexPoles));
        oldPhi = imag(complexPoles);
      % real ones
        oldA = sqrt(-oldA(imag(oldA) == 0));    

    epsilon = constrEps(G_BLA, G_var, s, oldA, oldAlpha, oldPhi, oldB);
    oldCost = epsilon' * epsilon;

%% Newton-Gauss method

    while 1
        
        [J, eps] = constrJacob(G_BLA, G_var, s, oldA, oldAlpha, oldPhi, oldB);

        J_re = [real(J); imag(J)];
        eps_re = [real(eps); imag(eps)];

        deltaTheta = -J_re\eps_re;

        thetaChange(itrCnt) = norm(deltaTheta);

        newA = oldA + deltaTheta(1:size(oldA, 1));
        newAlpha = oldAlpha + deltaTheta(size(oldA, 1)+1:size(oldA, 1)+size(oldAlpha, 1));
        newPhi = oldPhi + deltaTheta(size(oldA, 1)+size(oldAlpha, 1)+1:size(oldA, 1)+size(oldAlpha, 1)+size(oldPhi, 1));
        newB = oldB + deltaTheta(size(oldA, 1)+size(oldAlpha, 1)+size(oldPhi, 1)+1:end);

        epsilon = constrEps(G_BLA, G_var, s, newA, newAlpha, newPhi, newB);
    
        cost = epsilon' * epsilon;
        
        if ~isfinite(cost)
            error('ML iteration aborted: cost is NaN or Inf');
        end

        if cost >= oldCost
            % disp(['ML optimization finished after ', num2str(itrCnt), ' steps of Newton-Gauss method']);
            break;
        else
            oldCost = cost;
            oldA = newA;
            oldAlpha = newAlpha;
            oldPhi = newPhi;
            oldB = newB;
            itrCnt = itrCnt + 1;
        end

    end

    if show
        figure;
        plot(db(thetaChange), Color=[0.6350, 0.0780, 0.1840], LineWidth=2);
        xlabel('ML iteration');
        ylabel('Norm of \Delta\theta [dB]');
    end

    %% Levenberg-Marquardt

    [J, ~] = constrJacob(G_BLA, G_var, s, oldA, oldAlpha, oldPhi, oldB);

        J_re = [real(J); imag(J)];

    [~, S, ~] = svd(J_re);
    lambda = max(max(S))/100;

    itrCnt = 1;

    while 1

        [J, eps] = constrJacob(G_BLA, G_var, s, oldA, oldAlpha, oldPhi, oldB);

        J_re = [real(J); imag(J)];
        eps_re = [real(eps); imag(eps)];

        [U, S, V] = svd(J_re);
          % taking the rank deficiency of J into account
          S(size(S, 2), size(S, 2)) = 0;

        deltaTheta = -V * ((S'*S + lambda^2 * eye(size(S, 2)))\S') * U' * eps_re;

        % stop if relative parameter change gets low
        if (norm(deltaTheta)/(1+norm([oldA; oldAlpha; oldPhi; oldB]))) < 1e-12
            break;
        end

        newA = oldA + deltaTheta(1:size(oldA, 1));
        newAlpha = oldAlpha + deltaTheta(size(oldA, 1)+1:size(oldA, 1)+size(oldAlpha, 1));
        newPhi = oldPhi + deltaTheta(size(oldA, 1)+size(oldAlpha, 1)+1:size(oldA, 1)+size(oldAlpha, 1)+size(oldPhi, 1));
        newB = oldB + deltaTheta(size(oldA, 1)+size(oldAlpha, 1)+size(oldPhi, 1)+1:end);

        epsilon = constrEps(G_BLA, G_var, s, newA, newAlpha, newPhi, newB);

        cost = epsilon' * epsilon;

        if ~isfinite(cost)
            error('ML iteration aborted: cost is NaN or Inf');
        end

        if cost >= oldCost
            lambda = 10*lambda;
        else
            oldCost = cost;
            oldA = newA;
            oldAlpha = newAlpha;
            oldPhi = newPhi;
            oldB = newB;
            itrCnt = itrCnt + 1;
            lambda = 0.4*lambda;
        end

    end

    % Are contains real poles, Acplx complex ones
    Are = [];
    Acplx = [];
    if size(oldA, 1) ~= 0
        Are = -oldA.^2;
    end
    if size(oldAlpha, 1) ~= 0
        Acplx = [-oldAlpha.^2 + 1j*oldPhi; -oldAlpha.^2 - 1j*oldPhi];
    end
    A = [Are; Acplx];
    B = flipud(oldB);
    cost = oldCost;

    if ~isfinite(cost)
        error('Should never happen');
    end

%% Poles uncertainty computation
    thetaCov = inv(2*(J_re'*J_re));

    na = size(oldA, 1);
    nalpha = size(oldAlpha, 1);

    dp_dtheta = zeros(na+2*nalpha, na+2*nalpha);
    dp_dtheta(1:na, 1:na) = -2*diag(oldA);
    dp_dtheta(na+1:na+nalpha, na+1:na+nalpha) = -2*diag(oldAlpha);
    dp_dtheta(na+nalpha+1:na+2*nalpha, na+1:na+nalpha) = -2*diag(oldAlpha);
    dp_dtheta(na+1:na+nalpha, na+nalpha+1:na+2*nalpha) = 1j*eye(nalpha);
    dp_dtheta(na+nalpha+1:na+2*nalpha, na+nalpha+1:na+2*nalpha) = -1j*eye(nalpha);

      % split real and imaginary part
      dp_dtheta = [real(dp_dtheta); imag(dp_dtheta)];

    cov_p = dp_dtheta * thetaCov(1:na+2*nalpha, 1:na+2*nalpha) * dp_dtheta';

      % p_var 1 line for each pole
      % The first column contains the variance of the real part, the second
      % column contains the one the imaginary part.
    p_var = [diag(cov_p(1:Na, 1:Na)), diag(cov_p(Na+1:end, Na+1:end))];

end



function [jacob, eps] = constrJacob(G_BLA, G_var, s, oldA, oldAlpha, oldPhi, oldB)
    %% ----------- Jacobian matrix construction ---------------        
        % A:
            % Are is due to real poles, Acplx to complex ones
            Are = 1;
            Acplx = 1;
            if size(oldA, 1) ~= 0
                Are = prod(s + (oldA').^2, 2);
            end
            if size(oldAlpha, 1) ~= 0
                Acplx = prod(s.^2 + 2*s*(oldAlpha').^2 + (oldAlpha').^4 + (oldPhi').^2, 2);
            end
            A = Are .* Acplx;
        % B:
            B = s.^(0:size(oldB, 1)-1)*oldB;
        % e:
            e = A.*G_BLA-B;
        % d/dtheta (A, B):
            dA_dtheta_a = [];
            dA_dtheta_alpha = [];
            dA_dtheta_phi = [];
            if size(oldA, 1) ~= 0
                dA_dtheta_a = 2*A .* (oldA' ./ (s + (oldA').^2));
            end
            if size(oldAlpha, 1) ~= 0
                dA_dtheta_alpha = A .* ((2*s*oldAlpha' + 4*(oldAlpha').^3) ./ (s.^2 + 2*s*(oldAlpha').^2 + (oldAlpha').^4 + (oldPhi').^2));
                dA_dtheta_phi = A .* (2*oldPhi' ./ (s.^2 + 2*s*(oldAlpha').^2 + (oldAlpha').^4 + (oldPhi').^2));
            end
            dA_dtheta = [dA_dtheta_a, dA_dtheta_alpha, dA_dtheta_phi, zeros(size(G_BLA, 1), size(oldB, 1))];

            dB_dtheta = [zeros(size(G_BLA, 1), size(oldA, 1)+size(oldAlpha, 1)+size(oldPhi, 1)), s.^(0:size(oldB, 1)-1)];
        % d/dtheta (e):
            de_dthetaA = dA_dtheta(:, 1:size(oldA, 1)+size(oldAlpha, 1)+size(oldPhi, 1)).*G_BLA;
            de_dthetaB = -dB_dtheta(:, size(oldA, 1)+size(oldAlpha, 1)+size(oldPhi, 1)+1:end);
            de_dtheta = [de_dthetaA, de_dthetaB];
        % sigma_e:
            sigma_e = abs(A) .* sqrt(G_var);
        % d/dtheta (sigma_e^2):
            dsig_dtheta = 2*real(dA_dtheta.*conj(A)).*G_var;
    
        jacob = de_dtheta./sigma_e - e./(2*sigma_e.^3).*dsig_dtheta;

        eps = e./sigma_e;
end

function eps = constrEps(G_BLA, G_var, s, oldA, oldAlpha, oldPhi, oldB)
        % A:
            % Are is due to real poles, Acplx to complex ones
            Are = 1;
            Acplx = 1;
            if size(oldA, 1) ~= 0
                Are = prod(s + (oldA').^2, 2);
            end
            if size(oldAlpha, 1) ~= 0
                Acplx = prod(s.^2 + 2*s*(oldAlpha').^2 + (oldAlpha').^4 + (oldPhi').^2, 2);
            end
            A = Are .* Acplx;
        % B:
            B = s.^(0:size(oldB, 1)-1)*oldB;
        % e:
            e = A.*G_BLA-B;
        % sigma_e:
            sigma_e = abs(A) .* sqrt(G_var);

        eps = (e./sigma_e);
end



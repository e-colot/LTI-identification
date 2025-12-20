function [A, B, cost] = stableML(G_BLA, G_var, f, Na, Nb, show)
% [A, B, cost] = stableML(G_BLA, G_var, f, Na, Nb, [show])
% Determines the parameters A and B using a ML estimator that imposes 
% stable poles. The starting value for theta is computed using a GTLS 
% estimator
% 
% Na and Nb are the order of the denominator and numerator of the
% parametric transfer function
%
% show is 0 by default and shows how theta changes with the iterations
%
% Written by E. Colot on Dec 19 2025

    if nargin <= 5
        show = 0;
    end

    s = 1j*2*pi*f;
    [Agtls, Bgtls, ~] = GTLS(G_BLA, G_var, f, Na, Nb);
    oldA = roots(Agtls);
    oldB = flipud(Bgtls);
    oldCost = inf;
    itrCnt = 1;

    thetaChange = 0;

%% Newton-Gauss method

    while 1
        
        [J_re, eps_re] = constrJacob(G_BLA, G_var, s, oldA, oldB, Na, Nb);

        % removing low singular values
        tol = 1e-5;
        [U, S, V] = svd(J_re);
        S = S * diag((diag(S) > tol));
        J_re = U*S*V';

        % normalization
            S = diag(sqrt(sum(J_re.^2)));
            J_re = J_re * S;

        deltaTheta = -J_re\eps_re;

        % normalization
            deltaTheta = S * deltaTheta;

        thetaChange(itrCnt) = norm(deltaTheta);

        epsilon = constrEps(G_BLA, G_var, s, oldA + deltaTheta(1:Na), oldB + deltaTheta(Na+1:end), Nb);
    
        cost = epsilon' * epsilon;

        if cost >= oldCost
            % disp(['ML optimization finished after ', num2str(itrCnt), ' steps of Newton-Gauss method']);
            break;
        else
            oldCost = cost;
            oldA = oldA + deltaTheta(1:Na);
            oldB = oldB + deltaTheta(Na+1:end);
            itrCnt = itrCnt + 1;
        end

    end

    if show
        figure;
        plot(db(thetaChange), Color=[0.6350, 0.0780, 0.1840], LineWidth=2);
        xlabel('ML iteration');
        ylabel('Norm of \Delta\theta [dB]');
    end

    A = oldA;
    B = flipud(oldB);
    cost = oldCost;

end



function [jacob_re, eps_re] = constrJacob(G_BLA, G_var, s, oldA, oldB, Na, Nb)
    %% ----------- Jacobian matrix construction ---------------        
        % A:
            A = prod(s + (oldA.^2)', 2);
        % B:
            B = s.^(0:Nb)*oldB;
        % e:
            e = A.*G_BLA-B;
            e_re = [real(e); imag(e)];
        % d/dtheta (A):
            dA_dtheta = [2*A .* (oldA' ./ (s + (oldA.^2)')), zeros(size(G_BLA, 1), Nb+1)];
        % d/dtheta (e):
            de_r_dthetaAlpha = real(dA_dtheta(:, 1:Na)).*real(G_BLA) - imag(dA_dtheta(:, 1:Na)).*imag(G_BLA);
            de_i_dthetaAlpha = real(dA_dtheta(:, 1:Na)).*imag(G_BLA) + imag(dA_dtheta(:, 1:Na)).*real(G_BLA);
            de_r_dthetaB = -real(s.^(0:Nb));
            de_i_dthetaB = -imag(s.^(0:Nb));
            de_dtheta_re = [de_r_dthetaAlpha, de_r_dthetaB;
                            de_i_dthetaAlpha, de_i_dthetaB];
        % sigma_e:
            sigma_e = abs(A) .* sqrt(G_var);
        % d/dtheta (sigma_e^2):
            dsig_dtheta = 2*real(dA_dtheta.*conj(A)).*G_var;
    
        jacob_re = de_dtheta_re./repmat(sigma_e, 2, 1) - e_re./(2*repmat(sigma_e, 2, 1).^3).*repmat(dsig_dtheta, 2, 1);

        eps_re = (e_re./repmat(sigma_e, 2, 1));
end

function eps = constrEps(G_BLA, G_var, s, oldA, oldB, Nb)
        % A:
            A = prod(s + (oldA.^2)', 2);
        % B:
            B = s.^(0:Nb)*oldB;
        % e:
            e = A.*G_BLA-B;
        % sigma_e:
            sigma_e = abs(A) .* sqrt(G_var);

        eps = (e./sigma_e);
end


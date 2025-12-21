function [A, B, cost, K] = stableRealML(G_BLA, G_var, f, Na, Nb, show)
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
% Written by E. Colot on Dec 21 2025

    if nargin <= 5
        show = 0;
    end

    s = 1j*2*pi*f;
    [Agtls, Bgtls, ~] = GTLS(G_BLA, G_var, f, Na, Nb);
    oldA = roots(Agtls);
    oldB = flipud(Bgtls);
        % scale B to adapt
        oldB = oldB / Agtls(1);
    oldCost = inf;
    itrCnt = 1;

    thetaChange = 0;

      % impose negative roots
      oldA = -abs(real(oldA)) + 1j * imag(oldA);

    % split real and complex poles
    complexPoles = oldA(imag(oldA) > 0);
      % complex ones
        oldAlpha = abs(complexPoles);
        oldPhi = angle(complexPoles);
      % real ones
        oldA = oldA(imag(oldA) == 0);

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

        % check that the poles are still in the left half-plane
          kmax = 1;

            % check real poles
          for i = 1:size(oldA, 1)
              if newA(i) > 0
                  kmax = min(kmax, - (oldA(i)/(newA(i)-oldA(i))));
              end
          end
            % check complex poles
          for i = 1:size(oldAlpha, 1)
              oldPole = oldAlpha(i) * exp(1j * oldPhi(i));
              newPole = newAlpha(i) * exp(1j * newPhi(i));
              if real(newPole) > 0
                  kmax = min(kmax, - (real(oldPole)/(real(newPole-oldPole))));
              end
          end

        if kmax <= 0 || ~isfinite(kmax)
            error('Something went wrong...');
        end

        % safety margin
        kmax = 0.99*kmax;

        newA = oldA + kmax*(newA-oldA);
          % linearly changing amplitude and phase could result in a
          % positive pole -> going back to real-imag representation
        for i = 1:size(oldAlpha, 1)
              oldPole = oldAlpha(i) * exp(1j * oldPhi(i));
              newPole = newAlpha(i) * exp(1j * newPhi(i));
              newPole = oldPole + kmax*(newPole-oldPole);
              newAlpha(i) = abs(newPole);
              newPhi(i) = angle(newPole);
        end
        newB = oldB + kmax*(newB-oldB);

        % check for any overshoot in the right half-plane, meaning it is a
        % limit case and the optimization can be stopped
        stop = 0;
        if kmax < 1e-12
            % slow convergence, trying to push the poles in the right 
            % half-plane
            stop = 1;
        end
          for i = 1:size(oldA, 1)
              if newA(i) > 0
                  stop = 1;
              end
          end
            % check complex poles
          for i = 1:size(oldAlpha, 1)
              newPole = newAlpha(i) * exp(1j * newPhi(i));
              if real(newPole) > 0
                  stop = 1;
              end
          end
        if stop
            break;
        end

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

    % Are contains real poles, Acplx complex ones
    Are = [];
    Acplx = [];
    if size(oldA, 1) ~= 0
        Are = oldA;
    end
    if size(oldAlpha, 1) ~= 0
        Acplx = [oldAlpha .* exp(1j*oldPhi); oldAlpha .* exp(-1j*oldPhi)];
    end
    A = [Are; Acplx];
    B = flipud(oldB);
    cost = oldCost;

    if ~isfinite(cost)
        error('Should never happen');
    end

end



function [jacob, eps] = constrJacob(G_BLA, G_var, s, oldA, oldAlpha, oldPhi, oldB)
    %% ----------- Jacobian matrix construction ---------------        
        % A:
            % Are is due to real poles, Acplx to complex ones
            Are = 1;
            Acplx = 1;
            if size(oldA, 1) ~= 0
                Are = prod(s - oldA', 2);
            end
            if size(oldAlpha, 1) ~= 0
                Acplx = prod(s.^2 - 2*s*(oldAlpha' .* cos(oldPhi')) + (oldAlpha.^2)', 2);
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
                dA_dtheta_a = -A ./ (s - oldA');
            end
            if size(oldAlpha, 1) ~= 0
                dA_dtheta_alpha = A ./ (s.^2 - 2*s * oldAlpha' .* cos(oldPhi')  + (oldAlpha.^2)') .* (2*oldAlpha' - s*2.*cos(oldPhi'));
                dA_dtheta_phi = A ./ (s.^2 - 2*s * oldAlpha' .* cos(oldPhi')  + (oldAlpha.^2)') .* (s*2*oldAlpha'.*sin(oldPhi'));
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
                Are = prod(s - oldA', 2);
            end
            if size(oldAlpha, 1) ~= 0
                Acplx = prod(s.^2 - 2*s*(oldAlpha' .* cos(oldPhi')) + (oldAlpha.^2)', 2);
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



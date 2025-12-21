function [A, B, cost, K] = stableML(G_BLA, G_var, f, Na, Nb, show)
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
        % scale B to adapt
        oldB = oldB / Agtls(1);
    oldCost = inf;
    itrCnt = 1;

    thetaChange = 0;

      % impose negative roots
      oldA = -abs(real(oldA)) + 1j * imag(oldA);

%% Newton-Gauss method

    while 1
        
        [J, eps] = constrJacob(G_BLA, G_var, s, oldA, oldB, Na, Nb);

        deltaTheta = -J\eps;

        thetaChange(itrCnt) = norm(deltaTheta);

        newA = oldA + deltaTheta(1:Na);
        newB = oldB + [deltaTheta(Na+1:end)];

          % check that the poles are still in the left half-plane
          kmax = 1;
          for i = 1:Na
              if real(newA(i)) > 0
                  kmax = min(kmax, - (real(oldA(i))/(real(newA(i)-oldA(i)))));
              end
          end

        newA = oldA + kmax*(newA-oldA);
        newB = oldB + kmax*(newB-oldB);

        epsilon = constrEps(G_BLA, G_var, s, newA, newB, Nb);
    
        cost = epsilon' * epsilon;

        if cost >= oldCost
            % disp(['ML optimization finished after ', num2str(itrCnt), ' steps of Newton-Gauss method']);
            break;
        else
            oldCost = cost;
            oldA = newA;
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

    A = oldA;
    B = flipud(oldB);
    cost = oldCost;

end



function [jacob, eps] = constrJacob(G_BLA, G_var, s, oldA, oldB, Na, Nb)
    %% ----------- Jacobian matrix construction ---------------        
        % A:
            A = prod(s - oldA.', 2);
        % B:
            B = s.^(0:Nb)*oldB;
        % e:
            e = A.*G_BLA-B;
        % d/dtheta (A, B):
            dA_dtheta = [-A ./ (s - oldA.'), zeros(size(G_BLA, 1), Nb+1)];
        % d/dtheta (e):
            de_dthetaA = dA_dtheta(:, 1:Na).*G_BLA;
            de_dthetaB = -s.^(0:Nb);
            de_dtheta = [de_dthetaA, de_dthetaB];
        % sigma_e:
            sigma_e = abs(A) .* sqrt(G_var);
        % d/dtheta (sigma_e^2):
            dsig_dtheta = 2*real(dA_dtheta.*conj(A)).*G_var;
    
        jacob = de_dtheta./sigma_e - e./(2*sigma_e.^3).*dsig_dtheta;

        eps = e./sigma_e;
end

function eps = constrEps(G_BLA, G_var, s, oldA, oldB, Nb)
        % A:
            A = prod(s - oldA', 2);
        % B:
            B = s.^(0:Nb)*oldB;
        % e:
            e = A.*G_BLA-B;
        % sigma_e:
            sigma_e = abs(A) .* sqrt(G_var);

        eps = (e./sigma_e);
end



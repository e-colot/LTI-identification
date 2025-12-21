function [A, B, cost, thetaCov] = ML(G_BLA, G_var, f, Na, Nb, show)
% [A, B, cost, thetaCov] = ML(G_BLA, G_var, f, Na, Nb, [show])
% Determines the parameters A and B using a ML estimator. The starting
% value for theta is computed using a GTLS estimator
% 
% Na and Nb are the order of the denominator and numerator of the
% parametric transfer function
%
% show is 0 by default and shows how theta changes with the iterations
%
% Written by E. Colot on Dec 13 2025

    if nargin <= 5
        show = 0;
    end

    s = 1j*2*pi*f;
    [Agtls, Bgtls, ~] = GTLS(G_BLA, G_var, f, Na, Nb);
    theta0 = [flipud(Agtls); flipud(Bgtls)];
    
    oldTheta = theta0;
    oldCost = inf;
    itrCnt = 1;

    thetaChange = 0;

%% Newton-Gauss method

    while 1
        
        [jacob, epsilon] = constrJacob(G_BLA, G_var, s, oldTheta, Na, Nb);

        % For real theta:
            J_re = [real(jacob); imag(jacob)];
            eps_re = [real(epsilon); imag(epsilon)];

        % for conditioning, rms normalization of J (column-wise)
            
            % S = diag(1./sum(J_re.^2, 1));
            % J_re = J_re * S;

        deltaTheta = -J_re\eps_re;

        % deltaTheta = S * deltaTheta;

        thetaChange(itrCnt) = norm(deltaTheta);

        epsilon = constrEps(G_BLA, G_var, s, (oldTheta + deltaTheta)/norm(oldTheta + deltaTheta), Na, Nb);
    
        cost = epsilon' * epsilon;

        if cost >= oldCost
            % disp(['ML optimization finished after ', num2str(itrCnt), ' steps of Newton-Gauss method']);
            break;
        else
            oldCost = cost;
            oldTheta = oldTheta + deltaTheta;
            % to keep norm(theta) = 1
            oldTheta = oldTheta/norm(oldTheta);
            itrCnt = itrCnt + 1;
        end

    end

    if show
        figure;
        plot(db(thetaChange), Color=[0.6350, 0.0780, 0.1840], LineWidth=2);
        xlabel('ML iteration');
        ylabel('Norm of \Delta\theta [dB]');
    end

% %% Levenberg-Marquardt
% 
%     [jacob, ~] = constrJacob(G_BLA, G_var, s, oldTheta, Na, Nb);
%     [~, S, ~] = svd(jacob);
%     lambda = max(max(S))/100;
% 
%     oldCost = inf;
%     itrCnt = 0;
% 
%     while 1
% 
%         [jacob, epsilon] = constrJacob(G_BLA, G_var, s, oldTheta, Na, Nb);
% 
%         [U, S, V] = svd(jacob);
%           % taking the rank deficiency of J into account
%           S(size(S, 2), size(S, 2)) = 0;
% 
%         deltaTheta = -V * ((S.^2 + lambda^2 * [eye(size(S, 2)); zeros(size(S, 1)-size(S, 2), size(S, 2))])\S) * U.' * epsilon;
% 
%         epsilon = constrEps(G_BLA, G_var, s, (oldTheta + deltaTheta)/norm(oldTheta + deltaTheta), Na, Nb);
% 
%         cost = epsilon' * epsilon;
% 
%         if cost > oldCost
%             disp(['ML optimization finished after ', num2str(itrCnt), ' steps of Newton-Gauss method']);
%             break;
%         else
%             oldCost = cost;
%             oldTheta = oldTheta + deltaTheta;
%             % to keep norm(theta) = 1
%             oldTheta = oldTheta/norm(oldTheta);
%             itrCnt = itrCnt + 1;
%         end
% 
%     end


    A = flipud(oldTheta(1:Na+1));
    B = flipud(oldTheta(Na+2:end));
    cost = oldCost;

    thetaCov = inv(2*(J_re'*J_re));

end



function [jacob, eps] = constrJacob(G_BLA, G_var, s, oldTheta, Na, Nb)
    %% ----------- Jacobian matrix construction ---------------        
        % A:
            A = s.^(0:Na)*oldTheta(1:Na+1);
        % B:
            B = s.^(0:Nb)*oldTheta(Na+2:end);
        % e:
            e = A.*G_BLA-B;
        % d/dtheta (e):
            de_dtheta = s.^([(0:Na) (0:Nb)]);
            de_dtheta = de_dtheta.*[repmat(G_BLA, 1, Na+1), -ones(size(G_BLA, 1), Nb+1)];
        % sigma_e:
            sigma_e = abs(A) .* sqrt(G_var);
        % d/dtheta (sigma_e^2):
            dA_dtheta = [s.^(0:Na), zeros(size(G_BLA, 1), Nb+1)];
            dsig_dtheta = 2*real(dA_dtheta.*conj(A)).*G_var;
    
        jacob = de_dtheta./sigma_e - e./(2*sigma_e.^3).*dsig_dtheta;

        eps = (e./sigma_e);
end

function eps = constrEps(G_BLA, G_var, s, oldTheta, Na, Nb)
        % A:
            A = s.^(0:Na)*oldTheta(1:Na+1);
        % B:
            B = s.^(0:Nb)*oldTheta(Na+2:end);
        % e:
            e = A.*G_BLA-B;
        % sigma_e:
            sigma_e = abs(A) .* sqrt(G_var);

        eps = (e./sigma_e);
end


function [Na_opt, Nb_opt, A_opt, B_opt] = AIC(estimator, NaList, NbList, Ne)
% [Na_opt, Nb_opt, A_opt, B_opt] = AIC(estimator, NaList, NbList, Ne)
% finds the optimal order using the AIC criterion for the given estimator
% Ne is the number of points given to the estimator
%
% Written by E. Colot on Dec 12 2025

    cost_opt = inf;
    Na_opt = 0;
    Nb_opt = 0;
    A_opt  = [];
    B_opt  = [];

    % Storage for visualization
    AIC_map = NaN(length(NaList), length(NbList));

    for iNa = 1:length(NaList)
        Na = NaList(iNa);

        for iNb = 1:length(NbList)
            Nb = NbList(iNb);

            if Na < Nb
                % Stability requires Na >= Nb
                continue;
            end

            %disp(['Trying Na=', num2str(Na), ', Nb=', num2str(Nb)]);
            [A, B, cost] = estimator(Na, Nb);


            % AIC-like criterion
            AIC_cost = cost * (1 + (Na + Nb) / Ne);

            % Store
            AIC_map(iNa, iNb) = AIC_cost;

            if AIC_cost < cost_opt
                cost_opt = AIC_cost;
                Na_opt   = Na;
                Nb_opt   = Nb;
                A_opt    = A;
                B_opt    = B;
            end
        end
    end

    % ---- 3-D Visualization ----
    [Nb_grid, Na_grid] = meshgrid(NbList, NaList);
    AICdb = db(AIC_map);

    figure;
    surf(Nb_grid, Na_grid, AICdb, 'EdgeColor', 'none');
    xlabel('n_{b}');
    ylabel('n_{a}');
    zlabel('AIC cost [dB]');
    %title('AIC surface');
    colorbar;
    view(135, 30);
    grid on;
    hold on;

    % Mark optimum point
    plot3(Nb_opt, Na_opt, db(cost_opt), ...
          'ro', 'MarkerSize', 8, 'LineWidth', 2);

end

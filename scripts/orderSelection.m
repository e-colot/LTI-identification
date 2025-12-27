function [Na_opt, Nb_opt, A_opt, B_opt] = orderSelection(estimator, NaList, NbList, G_BLA, G_var, f)
% TODO
%
% Written by E. Colot on Dec 27 2025


%% Split training and validation sets

    setNumber = 5;

    cost_opt = inf;
    Na_opt = 0;
    Nb_opt = 0;

    % Storage for visualization
    residualCostMap = NaN(length(NaList), length(NbList));

    for iNa = 1:length(NaList)
        Na = NaList(iNa);

        for iNb = 1:length(NbList)
            Nb = NbList(iNb);

            if Na < Nb
                % Stability requires Na >= Nb
                break;
            end
            
            %disp(['Trying Na=', num2str(Na), ', Nb=', num2str(Nb)]);

            totalResidualCost = 0;

            for setCnt = 1:setNumber
                trainSet = setCnt:setNumber:length(f);
                validationSet = setdiff(1:length(f), trainSet);

                [~, ~, ~, costHandle] = estimator(G_BLA(trainSet), G_var(trainSet), f(trainSet), Na, Nb);
    
                totalResidualCost = totalResidualCost + costHandle(G_BLA(validationSet), G_var(validationSet), f(validationSet));
            end

            % Store
            residualCostMap(iNa, iNb) = totalResidualCost;

            if totalResidualCost < cost_opt
                cost_opt = totalResidualCost;
                Na_opt   = Na;
                Nb_opt   = Nb;
            end
        end
    end

    % ---- 3-D Visualization ----
    [Nb_grid, Na_grid] = meshgrid(NbList, NaList);

    figure;
    surf(Nb_grid, Na_grid, db(residualCostMap), 'EdgeColor', 'none');
    colormap(copper);
    xlabel('n_{b}');
    ylabel('n_{a}');
    zlabel('Residual Cost [dB]');
    %title('Residual Cost Surface');
    %colorbar;
    view(135, 30);
    grid on;
    hold on;

    % Mark optimum point
    plot3(Nb_opt, Na_opt, db(cost_opt), 'ro', 'MarkerSize', 8, 'LineWidth', 2);

%% Compute the optimal parameters on the whole dataset
    [A_opt, B_opt, ~, ~] = estimator(G_BLA, G_var, f, Na_opt, Nb_opt);


end



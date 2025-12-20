function [h1, h2] = plotTF(A, B, f, type)
% Plots G = B/A over the provided frequency vector
% returns a handle
%
% type can either be 
%  -'poly' if A = sum_k a_k s^k
%  -'mult' if A = prod_k (s+a_k^2)
%
% Written by E. Colot on Dec 11 2025

        s = 1j*2*pi*f;

        if strcmp(type, 'poly')
            A_eval = polyval(A, s);
            B_eval = polyval(B, s);
            
            G = B_eval./A_eval;
    
            subplot(211);
            h1 = plot(f, db(G), 'LineWidth', 1.5);
            subplot(212);
            h2 = plot(f, angle(G), 'LineWidth', 1.5);
        elseif strcmp(type, 'mult')
            A_eval = prod(s + (A.^2)', 2);
            B_eval = polyval(B, s);

            G = B_eval./A_eval;
    
            subplot(211);
            h1 = plot(f, db(G), 'LineWidth', 1.5);
            subplot(212);
            h2 = plot(f, angle(G), 'LineWidth', 1.5);
        else
            error('Wrong input argument "type"');
        end
end

function h = plotTF(A, B, f)
% Plots G = B/A over the provided frequency vector
% returns a handle
%
% Written by E. Colot on Dec 11 2025

        s = 1j*2*pi*f;
        A_eval = polyval(A, s);
        B_eval = polyval(B, s);
        
        G = B_eval./A_eval;

        h = plot(f, db(G), 'LineWidth', 1.5);

end

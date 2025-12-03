clear; close all; clc;

fs = 5e3;

  % number of poles
    n_a = 5;
  % number of zeros
    n_b = 4;

%% get a non parametric model

[f, G_ML, total_var] = robustMethod("robustMethod/full_5k", fs);

  % put f, G_ML and total_var as a column vectors
    if (size(f, 1) < size(f, 2))
        f = f';
    end
    if (size(G_ML, 1) < size(G_ML, 2))
        G_ML = G_ML';
    end
    if (size(total_var, 1) < size(total_var, 2))
        total_var = total_var';
    end
  % remove NaN's from G_ML (at Nyquist frequency);
    valid = ~isnan(G_ML);
    f       = f(valid);
    G_ML    = G_ML(valid);
    total_var = total_var(valid);

    

%% Initial estimate using GTLS

% J = [Y sY ... U sU ...]
% taking U = 1 and Y = G_ML

  % placing s in the matrix J
    J = repmat(1j*2*pi*f, 1, n_a+n_b+2);
    J = J.^([(0:n_a) (0:n_b)]);
  % placing G_ML in the matrix
    J = J.*[repmat(G_ML, 1, n_a+1) ones(size(G_ML, 1), n_b+1)];

% C_J = column covariance of J
    
  % placing sqrt(s) in the matrix J
    C_J = repmat(1j*2*pi*f, 1, n_a+1);
    C_J = C_J.^(0:n_a);
  % placing total_var in the matrix
    C_J = C_J .* repmat(total_var, 1, n_a+1);
  % the column with U have no variance (U = 1)
    C_J = [C_J zeros(size(G_ML, 1), n_b+1)];

% compute the gsvd

[U,V,X,C,S] = gsvd(J, sqrt(C_J));

[~, i] = min(diag(C)./diag(S));

Xinv = inv(X');
theta_GTLS = 1/S(i, i) * Xinv(:, i);





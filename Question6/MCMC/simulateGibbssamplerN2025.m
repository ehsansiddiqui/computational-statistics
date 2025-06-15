% simulateGibbssamplerN2025.m
% Modified to include MH in plots and loop over dimensions m = 2, 3, ..., 8

clear all;
close all;

% Set student number for reproducibility
studentnumber = 000611215;
randn('state', studentnumber);
rand('state', studentnumber);

% Parameters
n = 5000; % Number of samples
dims = 2:8; % Dimensions to test

for m = dims
    % Generate covariance matrix
    Z = randn(m, n);
    A = rand(m, m);
    Sigma = A * A';
    Kappa = inv(Sigma);
    
    %%%%%%%%%%%%%%%%%%%%% Transformation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y = A * Z;
    
    %%%%%%%%%%%%%%%%%%%%% Gibbs sampler %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X = zeros(m, n);
    B = cell(1, m);
    for k = 1:m
        notk = setdiff(1:m, k);
        B{k} = Sigma(k, notk) / Sigma(notk, notk);
    end
    Xi = zeros(m, 1);
    for i = 1:n
        for k = 1:m
            notk = setdiff(1:m, k);
            Xi(k) = B{k} * Xi(notk) + Z(k, i) / sqrt(Kappa(k, k));
        end
        X(1:m, i) = Xi;
    end
    
    %%%%%%%%%%%%%%%%%%%%% Metropolis-Hastings sampler %%%%%%%%%%%%%%%%%%%%%%
    V = zeros(m, n);
    x = V(1:m, 1);
    for i = 2:n
        y = NaN;
        while isnan(y)
            eta = rand(m, 1) * 2 - 1;
            y = x + eta;
            alfa = exp(x' * Kappa / 2 * x - y' * Kappa / 2 * y);
            U = rand;
            if alfa < U
                y = NaN;
            end
        end
        V(1:m, i) = y;
        x = y;
    end
    
    %%%%%%%%%%%%%%%%%%%%% Compute Sample Covariances %%%%%%%%%%%%%%%%%%%%%
    SigmaYhat = zeros(m, m);
    SigmaXhat = zeros(m, m);
    SigmaVhat = zeros(m, m);
    for i = 1:n
        SigmaYhat = SigmaYhat + Y(1:m, i) * Y(1:m, i)';
        SigmaXhat = SigmaXhat + X(1:m, i) * X(1:m, i)';
        SigmaVhat = SigmaVhat + V(1:m, i) * V(1:m, i)';
    end
    SigmaYhat = SigmaYhat / n;
    SigmaXhat = SigmaXhat / n;
    SigmaVhat = SigmaVhat / n;
    
    %%%%%%%%%%%%%%%%%%%%% Figure for m >= 2 %%%%%%%%%%%%%%%%%%%%%
    if m >= 2
        figure;
        hold on;
        plot(Y(1, :), Y(2, :), 'b.', 'DisplayName', 'Direct');
        plot(X(1, :), X(2, :), 'r.', 'DisplayName', 'Gibbs');
        plot(V(1, :), V(2, :), 'g.', 'DisplayName', 'MH');
        hold off;
        title(sprintf('Scatter Plot of First Two Dimensions, d = %d', m));
        xlabel('Dimension 1');
        ylabel('Dimension 2');
        legend('show');
        grid on;
    end
end
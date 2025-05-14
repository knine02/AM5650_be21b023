%% Generate Rössler data
tspan = 0:0.01:1000;
x0 = [1 1 1];
rossler = @(t, x)[-x(2)-x(3); x(1)+0.2*x(2); 0.2 + x(3)*(x(1) - 5.7)];
[~, X] = ode45(rossler, tspan, x0);
x = X(:,1);  % Choose x-component
x = x(1:10:end);  % Downsample

%% Parameters
m_list = 2:10;           % Embedding dimensions
tau = 10;               % Delay
epsilons = logspace(-5, 1, 50);  % Range of epsilon values
N = length(x);

%% Compute correlation sum C2 for each m and epsilon
log_eps = log(epsilons);
colors = lines(length(m_list));

figure; hold on;
for mi = 1:length(m_list)
    m = m_list(mi);
    % Reconstruct state space
    Y = zeros(N - (m - 1) * tau, m);
    for j = 1:m
        Y(:,j) = x((1:(N - (m - 1) * tau)) + (j - 1) * tau);
    end

    % Compute pairwise distances
    D = pdist2(Y, Y);
    D = D + eye(size(D)) * max(D(:));  % ignore self-distances

    % Correlation sum for each epsilon
    C2 = zeros(size(epsilons));
    for k = 1:length(epsilons)
        C2(k) = sum(D(:) < epsilons(k)) / (numel(D));
    end

    % Plot log-log
    plot(log_eps, log(C2), 'DisplayName', ['m = ' num2str(m)], ...
         'Color', colors(mi,:), 'LineWidth', 1.5);
end

xlabel('log(\epsilon)');
ylabel('log(C_2(\epsilon))');
legend('Location', 'northwest');
title('Correlation Sum vs. \epsilon for different m');
grid on;
%% Estimate D2 from linear region and plot D2 vs m
scaling_range = [0, 2];  % Choose log(epsilon) range for linear fit
D2_estimates = zeros(size(m_list));

for mi = 1:length(m_list)
    m = m_list(mi);

    % Reconstruct state space
    Y = zeros(N - (m - 1) * tau, m);
    for j = 1:m
        Y(:,j) = x((1:(N - (m - 1) * tau)) + (j - 1) * tau);
    end

    % Compute distances
    D = pdist2(Y, Y);
    D = D + eye(size(D)) * max(D(:));

    % Compute C2
    C2 = zeros(size(epsilons));
    for k = 1:length(epsilons)
        C2(k) = sum(D(:) < epsilons(k)) / (numel(D));
    end

    % Extract scaling region
    idx = log_eps > scaling_range(1) & log_eps < scaling_range(2);
    p = polyfit(log_eps(idx), log(C2(idx)), 1);
    D2_estimates(mi) = p(1);  % slope ≈ D2
end

%% Plot D2 vs m
figure;
plot(m_list, D2_estimates, '-o', 'LineWidth', 2);
xlabel('Embedding Dimension m');
ylabel('Estimated Correlation Dimension D_2');
title('Convergence of D_2 with Increasing m');
grid on;

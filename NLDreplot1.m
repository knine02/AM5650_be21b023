% Rössler system simulation and delay embedding visualization
clear; clc; close all;

% Parameters for Rössler system
a = 0.2; b = 0.2; c = 5.7;
dt = 0.05; T = 1000;
tspan = 0:dt:T;
x0 = [1, 1, 1];

% Define Rössler system
rossler = @(t, x)[-x(2) - x(3);
                  x(1) + a * x(2);
                  b + x(3) * (x(1) - c)];

% Integrate using ode45
[~, X] = ode45(rossler, tspan, x0);
x = X(:,1);  % Only use x-variable for time series

% Plot (a) Time series
figure('Position', [100, 100, 1200, 500]);
subplot(2,1,1)
plot(tspan, x, 'b')
xlabel('t')
ylabel('x')
title('(a)')
xlim([0, T])

% Embedding parameters
taus = [1, 10, 25, 35];  % Delay values
subplot_pos = 1;

% Plot (b) Delay embeddings
for i = 1:length(taus)
    tau = taus(i);
    % Create delay coordinates
    N = length(x) - 2 * tau;
    x1 = x(1:N);
    x2 = x((1:N) + tau);
    x3 = x((1:N) + 2 * tau);

    subplot(2, length(taus), length(taus) + i)
    plot3(x1, x2, x3, 'b', 'LineWidth', 0.5)
    view(2)
    axis tight
    axis square
    xlabel('x(t)')
    ylabel(['x(t + ', num2str(tau*dt), ')'])
    title(['\tau = ', num2str(tau)])
end

sgtitle('Rössler System and Delay Embeddings')

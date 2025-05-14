% Time vector
t = linspace(0, 20, 1000);  % time from 0 to 20 seconds

% Define time-varying frequency (nonstationarity)
f_t = 0.5 + 2 * (t / max(t));  % frequency increases from 0.5 to 2.5 Hz

% Compute phase by integrating the instantaneous frequency
phi_t = 2 * pi * cumtrapz(t, f_t);  % phase(t) = ∫ 2πf(t) dt

% Nonstationary sinusoidal signal
x = 20 * (sin(phi_t));

% Plot the signal
figure;
subplot(2,1,1)
plot(t, x, 'k')
xlabel('time')
ylabel('x')
title('Nonstationary Sinusoidal Signal')
xlim([0 20])

% Recurrence plot parameters
N = length(x);
eps = 1.5;  % distance threshold for recurrence

% Recurrence matrix
R = abs(x' - x) < eps;

% Plot the recurrence matrix
subplot(2,1,2)
imagesc(t, t, R)
colormap(flipud(gray))
axis square
xlabel('time')
ylabel('time')
title('Recurrence Plot')

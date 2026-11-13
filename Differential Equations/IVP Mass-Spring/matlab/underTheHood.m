%% underTheHood.m
% Two pictures of what a time stepper is doing. Run as is.
%   Part A. Euler and RK4 with the same fixed step, against the exact answer.
%   Part B. The step sizes ode45 actually chooses on a problem with a sudden kick.

clear; close all;

%% Part A. Euler vs. RK4, same step size
% Underdamped spring: y'' + 0.5 y' + y = 0, y(0) = 1, y'(0) = 0.
m = 1;  b = 0.5;  k = 1;  f = @(t) 0;
rhs = @(t, Y) [ Y(2); ( f(t) - b*Y(2) - k*Y(1) ) / m ];

% Exact solution of the underdamped free response
g  = b/(2*m);                 % decay rate
wd = sqrt(k/m - g^2);         % damped frequency
y0 = 1;  v0 = 0;
exact = @(t) exp(-g*t) .* ( y0*cos(wd*t) + (v0 + g*y0)/wd * sin(wd*t) );

h = 0.25;                     % one fixed step size for both methods
t = 0:h:20;
N = numel(t);

% Euler: follow the slope at the start of the step.
YE = zeros(2, N);  YE(:,1) = [y0; v0];
for n = 1:N-1
    YE(:,n+1) = YE(:,n) + h * rhs(t(n), YE(:,n));
end

% RK4: sample the slope four times across the step and average.
YR = zeros(2, N);  YR(:,1) = [y0; v0];
for n = 1:N-1
    k1 = rhs(t(n),       YR(:,n));
    k2 = rhs(t(n) + h/2, YR(:,n) + h/2*k1);
    k3 = rhs(t(n) + h/2, YR(:,n) + h/2*k2);
    k4 = rhs(t(n) + h,   YR(:,n) + h*k3);
    YR(:,n+1) = YR(:,n) + h/6 * (k1 + 2*k2 + 2*k3 + k4);
end

fprintf('Part A, step h = %g (%d steps)\n', h, N-1);
fprintf('   Euler max error: %.2e\n', max(abs(YE(1,:) - exact(t))));
fprintf('   RK4   max error: %.2e\n', max(abs(YR(1,:) - exact(t))));

tf = linspace(0, 20, 1000);
figure('Position', [100 100 1000 400]);
subplot(1, 2, 1);
plot(tf, exact(tf), 'k-', 'LineWidth', 1.5); hold on;
plot(t, YE(1,:), 'o-', 'MarkerSize', 4);
plot(t, YR(1,:), 's', 'MarkerSize', 4);
hold off; grid on;
xlabel('t'); ylabel('y(t)');
title(sprintf('Part A: same step h = %g', h));
legend('exact', 'Euler', 'RK4', 'Location', 'northeast');

%% Part B. Where ode45 puts its steps
% Lightly damped spring, released from y = 1 and then pushed at t = 10:
%   f(t) = 5 for 10 <= t <= 10.5, and 0 otherwise.
m = 1;  b = 0.2;  k = 1;
f = @(t) 5 * (t >= 10 & t <= 10.5);
rhs = @(t, Y) [ Y(2); ( f(t) - b*Y(2) - k*Y(1) ) / m ];

% Refine = 1 makes ode45 report only the steps it actually took.
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-9, 'Refine', 1);
[t, Y] = ode45(rhs, [0 30], [1; 0], opts);
hstep = diff(t);

fprintf('\nPart B: %d steps, smallest %.2e, largest %.2e\n', ...
        numel(hstep), min(hstep), max(hstep));

subplot(2, 2, 2);
plot(t, Y(:,1), '.-', 'MarkerSize', 6); grid on;
ylabel('y(t)');
title('Part B: a push at t = 10 (each dot is one step)');
subplot(2, 2, 4);
semilogy(t(1:end-1), hstep, '.-', 'MarkerSize', 6); grid on;
xlabel('t'); ylabel('step size h');

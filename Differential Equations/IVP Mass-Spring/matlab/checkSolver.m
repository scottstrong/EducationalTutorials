%% checkSolver.m
% Two ways to decide whether an ode45 answer can be trusted.
%   Check A. Compare against a case whose exact solution is known.
%   Check B. Change the tolerance and see whether the answer moves.
% Run as is. Nothing needs to be edited.

clear; close all;

%% Check A. Resonance, where the exact answer is known
% y'' + y = cos(t), y(0) = 0, y'(0) = 0 has exact solution y = (t/2) sin(t).
m = 1;  b = 0;  k = 1;  f = @(t) cos(t);
rhs   = @(t, Y) [ Y(2); ( f(t) - b*Y(2) - k*Y(1) ) / m ];
exact = @(t) t/2 .* sin(t);
tspan = [0 50];

tols   = [1e-3, 1e-6, 1e-9];      % 1e-3 is MATLAB's default RelTol
labels = cell(1, numel(tols));

figure('Position', [100 100 1000 400]);
subplot(1, 2, 1);
fprintf('Check A: resonance, compared with y = (t/2) sin t\n');
fprintf('   RelTol     output times    max |error|\n');
for j = 1:numel(tols)
    opts = odeset('RelTol', tols(j), 'AbsTol', tols(j)*1e-3);
    [t, Y] = ode45(rhs, tspan, [0; 0], opts);
    err = abs(Y(:,1) - exact(t));
    fprintf('   %7.0e   %10d       %10.2e\n', tols(j), numel(t), max(err));
    semilogy(t, err + eps, 'LineWidth', 1.2); hold on;
    labels{j} = sprintf('RelTol = %.0e', tols(j));
end
hold off; grid on;
xlabel('t'); ylabel('|ode45 - exact|');
ylim([1e-12 1]);
title('Check A: error against the exact solution');
legend(labels, 'Location', 'southeast');

%% Check B. Same problem, two tolerances, long run
% Undamped, unforced: y'' + y = 0, y(0) = 1, y'(0) = 0.
% The exact orbit is the unit circle, forever.
m = 1;  b = 0;  k = 1;  f = @(t) 0;
rhs   = @(t, Y) [ Y(2); ( f(t) - b*Y(2) - k*Y(1) ) / m ];
tspan = [0 1000];

[tL, YL] = ode45(rhs, tspan, [1; 0]);                                    % defaults
[tT, YT] = ode45(rhs, tspan, [1; 0], odeset('RelTol',1e-8,'AbsTol',1e-10)); % tight

fprintf('\nCheck B: undamped spring on [0, 1000], distance from origin at the end\n');
fprintf('   default tolerances: %.4f   (should be 1)\n', norm(YL(end,:)));
fprintf('   tight tolerances:   %.8f\n', norm(YT(end,:)));

subplot(1, 2, 2);
plot(tL, sqrt(YL(:,1).^2 + YL(:,2).^2), 'LineWidth', 1.5); hold on;
plot(tT, sqrt(YT(:,1).^2 + YT(:,2).^2), 'LineWidth', 1.5);
hold off; grid on; ylim([0.95 1.01]);
xlabel('t'); ylabel('distance from origin in the phase plane');
title('Check B: undamped spring, should stay at 1');
legend('default tolerances', 'RelTol = 1e-8', 'Location', 'southwest');

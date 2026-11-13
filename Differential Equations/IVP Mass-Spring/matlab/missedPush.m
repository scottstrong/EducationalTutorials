%% missedPush.m
% A short push applied to a spring at rest. With its default settings,
% ode45 can take a step so long that it jumps over the push and never sees it.
% Limiting the step size with MaxStep fixes it. Run as is.

clear; close all;

m = 1;  b = 0.2;  k = 1;
f = @(t) 5 * (t >= 10 & t <= 10.5);       % push of size 5 for half a second
rhs = @(t, Y) [ Y(2); ( f(t) - b*Y(2) - k*Y(1) ) / m ];
tspan = [0 30];
Y0 = [0; 0];                              % at rest, so nothing happens before t = 10

optsA = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
optsB = odeset('RelTol', 1e-6, 'AbsTol', 1e-9, 'MaxStep', 0.1);

[tA, YA] = ode45(rhs, tspan, Y0, optsA);
[tB, YB] = ode45(rhs, tspan, Y0, optsB);

fprintf('No MaxStep:      largest |y| = %.4f\n', max(abs(YA(:,1))));
fprintf('MaxStep = 0.1:   largest |y| = %.4f\n', max(abs(YB(:,1))));

figure('Position', [100 100 1000 400]);
plot(tA, YA(:,1), 'o-', 'LineWidth', 1.5, 'MarkerSize', 4); hold on;
plot(tB, YB(:,1), '-', 'LineWidth', 1.5);
hold off; grid on;
xlabel('t'); ylabel('y(t)');
title('A push at t = 10, starting from rest');
legend('default step control (missed the push)', 'MaxStep = 0.1', 'Location', 'northwest');

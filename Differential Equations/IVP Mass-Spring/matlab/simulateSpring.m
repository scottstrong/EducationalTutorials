%% simulateSpring.m
% Approximately solves the mass-spring initial value problem
%
%     m y'' + b y' + k y = f(t),     y(t0) = y0,   y'(t0) = v0
%
% with MATLAB's built-in adaptive solver ode45. Edit the values in
% Section 1, then press Run. Nothing below Section 1 needs to change.

clear; close all;

%% 1. Problem data (edit these)
m  = 1;           % mass
b  = 0.5;         % damping coefficient
k  = 1;           % spring constant
f  = @(t) 0;      % external forcing f(t); @(t) 0 means no forcing

y0 = 1;           % initial position  y(t0)
v0 = 0;           % initial velocity  y'(t0)
tspan = [0 30];   % [start time, end time]

%% 2. Rewrite as a first-order system
% The state is Y = [y; v] with v = y'. Then
%     y' = v
%     v' = ( f(t) - b*v - k*y ) / m
% Y(1) is position and Y(2) is velocity.
rhs = @(t, Y) [ Y(2);
                ( f(t) - b*Y(2) - k*Y(1) ) / m ];

%% 3. Solver settings
% RelTol and AbsTol set how much error ode45 allows in each step.
% (MATLAB's defaults are RelTol = 1e-3 and AbsTol = 1e-6, which are loose.)
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

%% 4. Solve
[t, Y] = ode45(rhs, tspan, [y0; v0], opts);
y = Y(:, 1);      % position at each output time
v = Y(:, 2);      % velocity at each output time

fprintf('ode45 returned %d output times on [%g, %g].\n', numel(t), tspan(1), tspan(end));
fprintf('Final state: y = %.6g, v = %.6g\n', y(end), v(end));

%% 5. Plot
figure('Position', [100 100 1000 400]);

subplot(1, 2, 1);
plot(t, y, 'LineWidth', 1.5);
grid on; xlim([tspan(1) tspan(end)]);
xlabel('t'); ylabel('y(t)');
title(sprintf('m = %g,  b = %g,  k = %g,  f = %s', m, b, k, func2str(f)));

subplot(1, 2, 2);
plot(y, v, 'LineWidth', 1.5); hold on;
plot(y(1), v(1), 'ko', 'MarkerFaceColor', 'k');   % starting point
hold off; grid on;
xlabel('y'); ylabel('v = y''');
title('Phase plane  (dot = start)');

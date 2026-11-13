# simulateSpring.R
# Approximately solves the mass-spring initial value problem
#
#     m y'' + b y' + k y = f(t),     y(t0) = y0,   y'(t0) = v0
#
# with the adaptive solver ode(..., method = "ode45") from the deSolve
# package (the same Dormand-Prince 4/5 pair as MATLAB's ode45). Edit the
# values in Section 1, then run. Nothing below Section 1 needs to change.
#
# R port of simulateSpring.m.  One-time setup: install.packages("deSolve")

library(deSolve)

## 1. Problem data (edit these)
m  <- 1                     # mass
b  <- 0.5                   # damping coefficient
k  <- 1                     # spring constant
f  <- function(t) 0         # external forcing f(t); function(t) 0 means no forcing

y0 <- 1                     # initial position  y(t0)
v0 <- 0                     # initial velocity  y'(t0)
tspan <- c(0, 30)           # c(start time, end time)

## 2. Rewrite as a first-order system
# The state is Y = c(y, v) with v = y'. Then
#     y' = v
#     v' = ( f(t) - b*v - k*y ) / m
# Y[1] is position and Y[2] is velocity. deSolve wants the derivative
# returned inside a list, and passes a third argument (parms) we do not use.
rhs <- function(t, Y, parms) {
  list(c(Y[2],
         (f(t) - b * Y[2] - k * Y[1]) / m))
}

## 3. Solver settings
# rtol and atol set how much error the solver allows in each step.
# (deSolve's defaults are rtol = 1e-6 and atol = 1e-6.)
# hmax is the largest step allowed. deSolve's default is the spacing of the
# output times, which would force tiny steps here; allow the whole interval
# and let the solver choose. (See missedPush.R for when to cap it.)
rtol <- 1e-8;  atol <- 1e-10;  hmax <- diff(tspan)

## 4. Solve
# deSolve reports the solution at the times you list; it still picks its own
# steps in between. Here: 1501 evenly spaced times across tspan.
times <- seq(tspan[1], tspan[2], length.out = 1501)
out <- ode(y = c(y = y0, v = v0), times = times, func = rhs, parms = NULL,
           method = "ode45", rtol = rtol, atol = atol, hmax = hmax)
t <- out[, "time"]
y <- out[, "y"]             # position at each output time
v <- out[, "v"]             # velocity at each output time

cat(sprintf("ode45 took %d steps on [%g, %g].\n",
            attr(out, "istate")[2], tspan[1], tspan[2]))
cat(sprintf("Final state: y = %.6g, v = %.6g\n", tail(y, 1), tail(v, 1)))

## 5. Plot
draw <- function() {
  par(mfrow = c(1, 2))
  plot(t, y, type = "l", lwd = 2, col = "blue", xlim = tspan,
       xlab = "t", ylab = "y(t)",
       main = sprintf("m = %g,  b = %g,  k = %g", m, b, k))
  grid()
  plot(y, v, type = "l", lwd = 2, col = "blue",
       xlab = "y", ylab = "v = y'", main = "Phase plane  (dot = start)")
  points(y[1], v[1], pch = 19)           # starting point
  grid()
}
if (interactive()) draw()
png("simulateSpring.png", width = 1000, height = 400); draw(); invisible(dev.off())

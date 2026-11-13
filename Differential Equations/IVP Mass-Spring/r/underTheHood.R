# underTheHood.R
# Two pictures of what a time stepper is doing. Run as is.
#   Part A. Euler and RK4 with the same fixed step, against the exact answer.
#   Part B. The step sizes ode45 actually chooses on a problem with a sudden kick.
#
# R port of underTheHood.m.  One-time setup: install.packages("deSolve")

library(deSolve)

## Part A. Euler vs. RK4, same step size
# Underdamped spring: y'' + 0.5 y' + y = 0, y(0) = 1, y'(0) = 0.
m <- 1;  b <- 0.5;  k <- 1;  f <- function(t) 0
rhs <- function(t, Y) c(Y[2], (f(t) - b * Y[2] - k * Y[1]) / m)

# Exact solution of the underdamped free response
g  <- b / (2 * m)                 # decay rate
wd <- sqrt(k / m - g^2)           # damped frequency
y0 <- 1;  v0 <- 0
exact <- function(t) exp(-g * t) * (y0 * cos(wd * t) + (v0 + g * y0) / wd * sin(wd * t))

h <- 0.25                         # one fixed step size for both methods
t <- seq(0, 20, by = h)
N <- length(t)

# Euler: follow the slope at the start of the step.
YE <- matrix(0, 2, N);  YE[, 1] <- c(y0, v0)
for (n in 1:(N - 1)) {
  YE[, n + 1] <- YE[, n] + h * rhs(t[n], YE[, n])
}

# RK4: sample the slope four times across the step and average.
YR <- matrix(0, 2, N);  YR[, 1] <- c(y0, v0)
for (n in 1:(N - 1)) {
  k1 <- rhs(t[n],         YR[, n])
  k2 <- rhs(t[n] + h / 2, YR[, n] + h / 2 * k1)
  k3 <- rhs(t[n] + h / 2, YR[, n] + h / 2 * k2)
  k4 <- rhs(t[n] + h,     YR[, n] + h * k3)
  YR[, n + 1] <- YR[, n] + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
}

cat(sprintf("Part A, step h = %g (%d steps)\n", h, N - 1))
cat(sprintf("   Euler max error: %.2e\n", max(abs(YE[1, ] - exact(t)))))
cat(sprintf("   RK4   max error: %.2e\n", max(abs(YR[1, ] - exact(t)))))

## Part B. Where ode45 puts its steps
# Lightly damped spring, released from y = 1 and then pushed at t = 10:
#   f(t) = 5 for 10 <= t <= 10.5, and 0 otherwise.
m <- 1;  b <- 0.2;  k <- 1
f <- function(t) 5 * (t >= 10 & t <= 10.5)
rhsB <- function(t, Y, parms) list(c(Y[2], (f(t) - b * Y[2] - k * Y[1]) / m))

# deSolve only reports the times you ask for, not the steps it took, so we
# log every time at which the right-hand side is evaluated and read the steps
# off that log. Each attempted step from t0 with size h evaluates the slope at
# t0 + h/5 and then t0 + 3h/10, so two consecutive calls a, b with the next
# call at a + 6(b - a) mark an attempt with h = 10(b - a), t0 = a - 2(b - a).
# An attempt was accepted when the next attempt starts at t0 + h.
calls <- numeric(0)
rhsLog <- function(t, Y, parms) { calls <<- c(calls, t); rhsB(t, Y, parms) }
invisible(ode(c(y = 1, v = 0), c(0, 30), rhsLog, NULL, method = "ode45",
              rtol = 1e-6, atol = 1e-9, hmax = 30))
d  <- diff(calls)
i  <- which(d[-length(d)] > 0 &
            abs(calls[-(1:2)] - calls[-(length(calls) - 0:1)] - 6 * d[-length(d)]) < 1e-9 * 30)
t0 <- calls[i] - 2 * d[i]                         # start of each attempt
ha <- 10 * d[i]                                   # size of each attempt
ok <- abs(c(t0[-1], 30) - (t0 + ha)) < 1e-9 * 30  # accepted attempts
tsteps <- c(t0[ok], 30)
hstep  <- diff(tsteps)
outB   <- ode(c(y = 1, v = 0), tsteps, rhsB, NULL, method = "ode45",
              rtol = 1e-6, atol = 1e-9, hmax = 30)

cat(sprintf("\nPart B: %d steps, smallest %.2e, largest %.2e\n",
            length(hstep), min(hstep), max(hstep)))

draw <- function() {
  layout(matrix(c(1, 1, 2, 3), 2, 2))
  tf <- seq(0, 20, length.out = 1000)
  plot(tf, exact(tf), type = "l", lwd = 2, ylim = range(YE[1, ]), xlab = "t", ylab = "y(t)",
       main = sprintf("Part A: same step h = %g", h))
  lines(t, YE[1, ], type = "o", pch = 1, cex = 0.6, col = "blue")
  points(t, YR[1, ], pch = 0, cex = 0.6, col = "red")
  grid()
  legend("topright", c("exact", "Euler", "RK4"), col = c("black", "blue", "red"),
         lty = c(1, 1, NA), pch = c(NA, 1, 0), bty = "n")

  plot(tsteps, outB[, "y"], type = "o", pch = 20, cex = 0.5, xlab = "", ylab = "y(t)",
       main = "Part B: a push at t = 10 (each dot is one step)")
  grid()
  plot(head(tsteps, -1), hstep, type = "o", pch = 20, cex = 0.5, log = "y",
       xlab = "t", ylab = "step size h")
  grid()
}
if (interactive()) draw()
png("underTheHood.png", width = 1000, height = 400); draw(); invisible(dev.off())

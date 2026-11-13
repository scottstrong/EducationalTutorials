# checkSolver.R
# Two ways to decide whether an ode45 answer can be trusted.
#   Check A. Compare against a case whose exact solution is known.
#   Check B. Change the tolerance and see whether the answer moves.
# Run as is. Nothing needs to be edited.
#
# R port of checkSolver.m.  One-time setup: install.packages("deSolve")

library(deSolve)

spring <- function(m, b, k, f) {
  # Right-hand side of the first-order system for m y'' + b y' + k y = f(t)
  function(t, Y, parms) list(c(Y[2], (f(t) - b * Y[2] - k * Y[1]) / m))
}

## Check A. Resonance, where the exact answer is known
# y'' + y = cos(t), y(0) = 0, y'(0) = 0 has exact solution y = (t/2) sin(t).
rhs   <- spring(1, 0, 1, cos)
exact <- function(t) t / 2 * sin(t)
times <- seq(0, 50, by = 0.05)
# hmax = the whole interval: otherwise deSolve caps the step at the output
# spacing, and every tolerance would give the same (tiny-step) answer.

tols <- c(1e-3, 1e-6, 1e-9)       # 1e-3 is MATLAB's default RelTol
errs <- list()
cat("Check A: resonance, compared with y = (t/2) sin t\n")
cat("   rtol      steps taken    max |error|\n")
for (j in seq_along(tols)) {
  out <- ode(c(y = 0, v = 0), times, rhs, NULL, method = "ode45",
             rtol = tols[j], atol = tols[j] * 1e-3, hmax = 50)
  errs[[j]] <- abs(out[, "y"] - exact(times))
  cat(sprintf("   %7.0e   %10d       %10.2e\n",
              tols[j], attr(out, "istate")[2], max(errs[[j]])))
}

## Check B. Same problem, two tolerances, long run
# Undamped, unforced: y'' + y = 0, y(0) = 1, y'(0) = 0.
# The exact orbit is the unit circle, forever.
rhs    <- spring(1, 0, 1, function(t) 0)
timesB <- seq(0, 1000, by = 0.5)

# "Loose" = MATLAB's defaults, so the comparison matches the tutorial.
outL <- ode(c(y = 1, v = 0), timesB, rhs, NULL, method = "ode45", rtol = 1e-3, atol = 1e-6, hmax = 1000)
outT <- ode(c(y = 1, v = 0), timesB, rhs, NULL, method = "ode45", rtol = 1e-8, atol = 1e-10, hmax = 1000)
rL <- sqrt(outL[, "y"]^2 + outL[, "v"]^2)
rT <- sqrt(outT[, "y"]^2 + outT[, "v"]^2)

cat("\nCheck B: undamped spring on [0, 1000], distance from origin at the end\n")
cat(sprintf("   loose tolerances:   %.4f   (should be 1)\n", tail(rL, 1)))
cat(sprintf("   tight tolerances:   %.8f\n", tail(rT, 1)))

draw <- function() {
  par(mfrow = c(1, 2))
  cols <- c("blue", "red", "darkgreen")
  plot(NULL, xlim = c(0, 50), ylim = c(1e-12, 1), log = "y",
       xlab = "t", ylab = "|ode45 - exact|",
       main = "Check A: error against the exact solution")
  for (j in seq_along(tols)) lines(times, errs[[j]] + .Machine$double.eps, col = cols[j])
  grid()
  legend("bottomright", sprintf("rtol = %.0e", tols), col = cols, lty = 1, bty = "n")

  plot(timesB, rL, type = "l", lwd = 2, col = "blue", ylim = c(0.95, 1.01),
       xlab = "t", ylab = "distance from origin in the phase plane",
       main = "Check B: undamped spring, should stay at 1")
  lines(timesB, rT, lwd = 2, col = "red")
  grid()
  legend("bottomleft", c("rtol = 1e-3 (MATLAB default)", "rtol = 1e-8"),
         col = c("blue", "red"), lwd = 2, bty = "n")
}
if (interactive()) draw()
png("checkSolver.png", width = 1000, height = 400); draw(); invisible(dev.off())

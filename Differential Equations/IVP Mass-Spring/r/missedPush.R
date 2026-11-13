# missedPush.R
# A short push applied to a spring at rest. With a loose limit on the step
# size, ode45 can take a step so long that it jumps over the push and never
# sees it. Limiting the step size with hmax fixes it. Run as is.
#
# R port of missedPush.m.  One-time setup: install.packages("deSolve")
#
# Note: deSolve's default hmax is the largest gap in `times`. Asking for
# output on a fine grid therefore already protects you; asking only for the
# two endpoints, as below, lets a single step span the whole interval.

library(deSolve)

m <- 1;  b <- 0.2;  k <- 1
f <- function(t) 5 * (t >= 10 & t <= 10.5)     # push of size 5 for half a second
rhs <- function(t, Y, parms) list(c(Y[2], (f(t) - b * Y[2] - k * Y[1]) / m))
Y0 <- c(y = 0, v = 0)                          # at rest, so nothing happens before t = 10

outA <- ode(Y0, c(0, 30), rhs, NULL, method = "ode45", rtol = 1e-6, atol = 1e-9)
outB <- ode(Y0, c(0, 30), rhs, NULL, method = "ode45", rtol = 1e-6, atol = 1e-9, hmax = 0.1)

# To draw the curves, re-run on a grid but keep each run's step limit.
tt <- seq(0, 30, by = 0.05)
curveA <- ode(Y0, tt, rhs, NULL, method = "ode45", rtol = 1e-6, atol = 1e-9, hmax = 30)
curveB <- ode(Y0, tt, rhs, NULL, method = "ode45", rtol = 1e-6, atol = 1e-9, hmax = 0.1)

cat(sprintf("No hmax:         largest |y| = %.4f\n", max(abs(curveA[, "y"]))))
cat(sprintf("hmax = 0.1:      largest |y| = %.4f\n", max(abs(curveB[, "y"]))))

draw <- function() {
  plot(tt, curveA[, "y"], type = "l", lwd = 2, col = "blue",
       ylim = range(curveB[, "y"]), xlab = "t", ylab = "y(t)",
       main = "A push at t = 10, starting from rest")
  lines(tt, curveB[, "y"], lwd = 2, col = "red")
  grid()
  legend("topleft", c("no step limit (missed the push)", "hmax = 0.1"),
         col = c("blue", "red"), lwd = 2, bty = "n")
}
if (interactive()) draw()
png("missedPush.png", width = 1000, height = 400); draw(); invisible(dev.off())

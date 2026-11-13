"""underTheHood.py
Two pictures of what a time stepper is doing. Run as is.
  Part A. Euler and RK4 with the same fixed step, against the exact answer.
  Part B. The step sizes solve_ivp actually chooses on a problem with a sudden kick.

Python port of underTheHood.m.  Needs numpy, scipy, matplotlib.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# %% Part A. Euler vs. RK4, same step size
# Underdamped spring: y'' + 0.5 y' + y = 0, y(0) = 1, y'(0) = 0.
m, b, k = 1.0, 0.5, 1.0
f = lambda t: 0.0
rhs = lambda t, Y: np.array([Y[1], (f(t) - b * Y[1] - k * Y[0]) / m])

# Exact solution of the underdamped free response
g = b / (2 * m)                    # decay rate
wd = np.sqrt(k / m - g**2)         # damped frequency
y0, v0 = 1.0, 0.0
exact = lambda t: np.exp(-g * t) * (y0 * np.cos(wd * t) + (v0 + g * y0) / wd * np.sin(wd * t))

h = 0.25                           # one fixed step size for both methods
t = np.arange(0, 20 + h / 2, h)
N = len(t)

# Euler: follow the slope at the start of the step.
YE = np.zeros((2, N)); YE[:, 0] = [y0, v0]
for n in range(N - 1):
    YE[:, n + 1] = YE[:, n] + h * rhs(t[n], YE[:, n])

# RK4: sample the slope four times across the step and average.
YR = np.zeros((2, N)); YR[:, 0] = [y0, v0]
for n in range(N - 1):
    k1 = rhs(t[n],         YR[:, n])
    k2 = rhs(t[n] + h / 2, YR[:, n] + h / 2 * k1)
    k3 = rhs(t[n] + h / 2, YR[:, n] + h / 2 * k2)
    k4 = rhs(t[n] + h,     YR[:, n] + h * k3)
    YR[:, n + 1] = YR[:, n] + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

print(f"Part A, step h = {h:g} ({N - 1} steps)")
print(f"   Euler max error: {np.max(np.abs(YE[0] - exact(t))):.2e}")
print(f"   RK4   max error: {np.max(np.abs(YR[0] - exact(t))):.2e}")

fig = plt.figure(figsize=(10, 4))
ax = fig.add_subplot(1, 2, 1)
tf = np.linspace(0, 20, 1000)
ax.plot(tf, exact(tf), "k-", lw=1.5, label="exact")
ax.plot(t, YE[0], "o-", ms=4, label="Euler")
ax.plot(t, YR[0], "s", ms=4, label="RK4")
ax.grid(True); ax.set_xlabel("t"); ax.set_ylabel("y(t)")
ax.set_title(f"Part A: same step h = {h:g}")
ax.legend(loc="upper right")

# %% Part B. Where solve_ivp puts its steps
# Lightly damped spring, released from y = 1 and then pushed at t = 10:
#   f(t) = 5 for 10 <= t <= 10.5, and 0 otherwise.
m, b, k = 1.0, 0.2, 1.0
f = lambda t: 5.0 * (10 <= t <= 10.5)
rhs = lambda t, Y: [Y[1], (f(t) - b * Y[1] - k * Y[0]) / m]

# Without t_eval, sol.t holds exactly the steps the solver took.
sol = solve_ivp(rhs, (0, 30), [1, 0], method="RK45", rtol=1e-6, atol=1e-9)
hstep = np.diff(sol.t)

print(f"\nPart B: {len(hstep)} steps, smallest {hstep.min():.2e}, largest {hstep.max():.2e}")

ax = fig.add_subplot(2, 2, 2)
ax.plot(sol.t, sol.y[0], ".-", ms=4); ax.grid(True)
ax.set_ylabel("y(t)")
ax.set_title("Part B: a push at t = 10 (each dot is one step)")
ax = fig.add_subplot(2, 2, 4)
ax.semilogy(sol.t[:-1], hstep, ".-", ms=4); ax.grid(True)
ax.set_xlabel("t"); ax.set_ylabel("step size h")

fig.tight_layout()
fig.savefig("underTheHood.png", dpi=150)
plt.show()

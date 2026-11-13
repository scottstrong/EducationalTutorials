"""missedPush.py
A short push applied to a spring at rest. With its default settings,
solve_ivp can take a step so long that it jumps over the push and never sees it.
Limiting the step size with max_step fixes it. Run as is.

Python port of missedPush.m.  Needs numpy, scipy, matplotlib.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

m, b, k = 1.0, 0.2, 1.0
f = lambda t: 5.0 * (10 <= t <= 10.5)       # push of size 5 for half a second
rhs = lambda t, Y: [Y[1], (f(t) - b * Y[1] - k * Y[0]) / m]
tspan = (0, 30)
Y0 = [0, 0]                                 # at rest, so nothing happens before t = 10

# SciPy's default max_step is infinity, so nothing stops a huge step.
solA = solve_ivp(rhs, tspan, Y0, method="RK45", rtol=1e-6, atol=1e-9)
solB = solve_ivp(rhs, tspan, Y0, method="RK45", rtol=1e-6, atol=1e-9, max_step=0.1)

print(f"No max_step:      largest |y| = {np.max(np.abs(solA.y[0])):.4f}")
print(f"max_step = 0.1:   largest |y| = {np.max(np.abs(solB.y[0])):.4f}")

plt.figure(figsize=(10, 4))
plt.plot(solA.t, solA.y[0], "o-", lw=1.5, ms=4, label="default step control (missed the push)")
plt.plot(solB.t, solB.y[0], "-", lw=1.5, label="max_step = 0.1")
plt.grid(True); plt.xlabel("t"); plt.ylabel("y(t)")
plt.title("A push at t = 10, starting from rest")
plt.legend(loc="upper left")
plt.tight_layout()
plt.savefig("missedPush.png", dpi=150)
plt.show()

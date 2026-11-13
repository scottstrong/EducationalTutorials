"""checkSolver.py
Two ways to decide whether a solve_ivp answer can be trusted.
  Check A. Compare against a case whose exact solution is known.
  Check B. Change the tolerance and see whether the answer moves.
Run as is. Nothing needs to be edited.

Python port of checkSolver.m.  Needs numpy, scipy, matplotlib.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def spring(m, b, k, f):
    """Right-hand side of the first-order system for m y'' + b y' + k y = f(t)."""
    return lambda t, Y: [Y[1], (f(t) - b * Y[1] - k * Y[0]) / m]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

# %% Check A. Resonance, where the exact answer is known
# y'' + y = cos(t), y(0) = 0, y'(0) = 0 has exact solution y = (t/2) sin(t).
rhs = spring(1, 0, 1, np.cos)
exact = lambda t: t / 2 * np.sin(t)
tspan = (0, 50)

tols = [1e-3, 1e-6, 1e-9]          # 1e-3 is SciPy's (and MATLAB's) default rtol
print("Check A: resonance, compared with y = (t/2) sin t")
print("   rtol      steps taken    max |error|")
for tol in tols:
    sol = solve_ivp(rhs, tspan, [0, 0], method="RK45", rtol=tol, atol=tol * 1e-3)
    err = np.abs(sol.y[0] - exact(sol.t))
    print(f"   {tol:7.0e}   {len(sol.t) - 1:10d}       {err.max():10.2e}")
    ax1.semilogy(sol.t, err + np.finfo(float).eps, lw=1.2, label=f"rtol = {tol:.0e}")
ax1.grid(True); ax1.set_ylim(1e-12, 1)
ax1.set_xlabel("t"); ax1.set_ylabel("|solve_ivp - exact|")
ax1.set_title("Check A: error against the exact solution")
ax1.legend(loc="lower right")

# %% Check B. Same problem, two tolerances, long run
# Undamped, unforced: y'' + y = 0, y(0) = 1, y'(0) = 0.
# The exact orbit is the unit circle, forever.
rhs = spring(1, 0, 1, lambda t: 0.0)
tspan = (0, 1000)

solL = solve_ivp(rhs, tspan, [1, 0], method="RK45")                          # defaults
solT = solve_ivp(rhs, tspan, [1, 0], method="RK45", rtol=1e-8, atol=1e-10)    # tight

print("\nCheck B: undamped spring on [0, 1000], distance from origin at the end")
print(f"   default tolerances: {np.linalg.norm(solL.y[:, -1]):.4f}   (should be 1)")
print(f"   tight tolerances:   {np.linalg.norm(solT.y[:, -1]):.8f}")

ax2.plot(solL.t, np.hypot(*solL.y), lw=1.5, label="default tolerances")
ax2.plot(solT.t, np.hypot(*solT.y), lw=1.5, label="rtol = 1e-8")
ax2.grid(True); ax2.set_ylim(0.95, 1.01)
ax2.set_xlabel("t"); ax2.set_ylabel("distance from origin in the phase plane")
ax2.set_title("Check B: undamped spring, should stay at 1")
ax2.legend(loc="lower left")

fig.tight_layout()
fig.savefig("checkSolver.png", dpi=150)
plt.show()

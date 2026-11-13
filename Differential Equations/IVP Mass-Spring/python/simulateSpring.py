"""simulateSpring.py
Approximately solves the mass-spring initial value problem

    m y'' + b y' + k y = f(t),     y(t0) = y0,   y'(t0) = v0

with SciPy's built-in adaptive solver solve_ivp (method "RK45", the same
Dormand-Prince 4/5 pair as MATLAB's ode45). Edit the values in Section 1,
then run. Nothing below Section 1 needs to change.

Python port of simulateSpring.m.  Needs numpy, scipy, matplotlib.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# %% 1. Problem data (edit these)
m = 1.0            # mass
b = 0.5            # damping coefficient
k = 1.0            # spring constant
f = lambda t: 0.0  # external forcing f(t); lambda t: 0.0 means no forcing

y0 = 1.0           # initial position  y(t0)
v0 = 0.0           # initial velocity  y'(t0)
tspan = (0, 30)    # (start time, end time)

# %% 2. Rewrite as a first-order system
# The state is Y = [y, v] with v = y'. Then
#     y' = v
#     v' = ( f(t) - b*v - k*y ) / m
# Y[0] is position and Y[1] is velocity (Python counts from 0).
def rhs(t, Y):
    return [Y[1],
            (f(t) - b * Y[1] - k * Y[0]) / m]

# %% 3. Solver settings
# rtol and atol set how much error the solver allows in each step.
# (SciPy's defaults are rtol = 1e-3 and atol = 1e-6, which are loose.)
opts = dict(method="RK45", rtol=1e-8, atol=1e-10)

# %% 4. Solve
# dense_output + t_eval gives a smooth curve to plot; the solver still
# chooses its own steps (sol.t would hold only the steps it took without t_eval).
t_eval = np.linspace(tspan[0], tspan[1], 1501)
sol = solve_ivp(rhs, tspan, [y0, v0], t_eval=t_eval, **opts)
t = sol.t
y = sol.y[0]       # position at each output time
v = sol.y[1]       # velocity at each output time

print(f"solve_ivp took {sol.nfev} slope evaluations on [{tspan[0]:g}, {tspan[1]:g}].")
print(f"Final state: y = {y[-1]:.6g}, v = {v[-1]:.6g}")

# %% 5. Plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

ax1.plot(t, y, lw=1.5)
ax1.grid(True); ax1.set_xlim(tspan)
ax1.set_xlabel("t"); ax1.set_ylabel("y(t)")
ax1.set_title(f"m = {m:g},  b = {b:g},  k = {k:g}")

ax2.plot(y, v, lw=1.5)
ax2.plot(y[0], v[0], "ko")                 # starting point
ax2.grid(True)
ax2.set_xlabel("y"); ax2.set_ylabel("v = y'")
ax2.set_title("Phase plane  (dot = start)")

fig.tight_layout()
fig.savefig("simulateSpring.png", dpi=150)
plt.show()

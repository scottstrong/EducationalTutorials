(* simulateSpring.wl
   Approximately solves the mass-spring initial value problem

       m y'' + b y' + k y = f(t),     y(t0) = y0,   y'(t0) = v0

   with Mathematica's NDSolve, told to use an adaptive explicit Runge-Kutta
   4/5 pair (the same kind of method as MATLAB's ode45). Edit the values in
   Section 1, then evaluate the whole file. Nothing below Section 1 needs to change.

   Wolfram port of simulateSpring.m.  Numerical only. *)

Quiet@Check[SetDirectory[NotebookDirectory[]], Null];   (* Export writes to Directory[] *)

(* ::Section:: *)
(* 1. Problem data (edit these) *)
m  = 1.;            (* mass *)
b  = 0.5;           (* damping coefficient *)
k  = 1.;            (* spring constant *)
f[t_?NumericQ] := 0.;  (* external forcing f(t); 0. means no forcing *)
(* The ?NumericQ makes NDSolve treat f as a black box it can only sample,
   exactly like ode45 does, instead of reading the formula. *)

y0 = 1.;            (* initial position  y(t0) *)
v0 = 0.;            (* initial velocity  y'(t0) *)
tspan = {0., 30.};  (* {start time, end time} *)

(* ::Section:: *)
(* 2. Rewrite as a first-order system *)
(* The state is {y, v} with v = y'. Then
       y' = v
       v' = ( f(t) - b v - k y ) / m                                       *)
rhs[t_, {yy_, vv_}] := {vv, (f[t] - b vv - k yy)/m};

(* ::Section:: *)
(* 3. Solver settings *)
(* PrecisionGoal -> 8 is RelTol = 1e-8; AccuracyGoal -> 10 is AbsTol = 1e-10. *)
opts = {Method -> {"ExplicitRungeKutta", "DifferenceOrder" -> 5},
        PrecisionGoal -> 8, AccuracyGoal -> 10};

(* ::Section:: *)
(* 4. Solve *)
Ysol = Y /. First@NDSolve[
   {Y'[t] == rhs[t, Y[t]], Y[tspan[[1]]] == {y0, v0}},
   Y, {t, tspan[[1]], tspan[[2]]}, Sequence @@ opts];
y[t_?NumericQ] := Ysol[t][[1]];   (* position at time t *)
v[t_?NumericQ] := Ysol[t][[2]];   (* velocity at time t *)

nSteps = Length[Ysol["Grid"]] - 1;   (* Grid = the times NDSolve stepped to *)
Print["NDSolve took ", nSteps, " steps on [", tspan[[1]], ", ", tspan[[2]], "]."];
Print["Final state: y = ", ToString[y[tspan[[2]]]], ", v = ", ToString[v[tspan[[2]]]]];

(* ::Section:: *)
(* 5. Plot *)
plot = GraphicsRow[{
   Plot[y[t], {t, tspan[[1]], tspan[[2]]}, PlotRange -> All,
     GridLines -> Automatic, Frame -> True, FrameLabel -> {"t", "y(t)"},
     PlotLabel -> StringForm["m = ``,  b = ``,  k = ``", m, b, k]],
   Show[
     ParametricPlot[{y[t], v[t]}, {t, tspan[[1]], tspan[[2]]}, PlotRange -> All,
       AspectRatio -> 1/GoldenRatio, GridLines -> Automatic, Frame -> True,
       FrameLabel -> {"y", "v = y'"}, PlotLabel -> "Phase plane  (dot = start)"],
     Graphics[{Black, PointSize[0.03], Point[{y0, v0}]}]]   (* starting point *)
   }, ImageSize -> 1000];

Export["simulateSpring.png", plot];
plot                                     (* evaluate on its own line to display *)

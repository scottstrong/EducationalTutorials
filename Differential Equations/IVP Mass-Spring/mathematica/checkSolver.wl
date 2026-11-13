(* checkSolver.wl
   Two ways to decide whether an NDSolve answer can be trusted.
     Check A. Compare against a case whose exact solution is known.
     Check B. Change the tolerance and see whether the answer moves.
   Run as is. Nothing needs to be edited.

   Wolfram port of checkSolver.m.  Numerical only.
   Tolerances: PrecisionGoal -> p means RelTol = 10^-p, AccuracyGoal -> a means AbsTol = 10^-a. *)

Quiet@Check[SetDirectory[NotebookDirectory[]], Null];   (* Export writes to Directory[] *)

(* Batch-safe scientific notation, e.g. sci[0.0000521] -> "5.21e-05" *)
sci[x_?NumericQ, d_: 2] := If[x == 0, "0", Module[{e = Floor[Log10[Abs[x]]], mm},
   mm = x/10.^e;
   If[d == 0, ToString[Round[mm]], ToString@NumberForm[mm, {d + 1, d}]] <>
    "e" <> If[e < 0, "-", "+"] <> IntegerString[Abs[e], 10, 2]]];

(* Solve m y'' + b y' + k y = f(t) as a first-order system; return the solution
   (a vector-valued interpolating function of t). f must be a black-box function. *)
(* The {yy_, vv_} pattern only matches a numeric state, so NDSolve sees the
   system as a black box it can sample, like ode45. *)
springRHS[m_, b_, k_, f_][t_, {yy_, vv_}] := {vv, (f[t] - b vv - k yy)/m};
solveSpring[m_, b_, k_, f_, Y0_, {t0_, t1_}, relTol_, absTol_] := Module[{Y, t},
  Y /. First@NDSolve[
     {Y'[t] == springRHS[m, b, k, f][t, Y[t]], Y[t0] == Y0},
     Y, {t, t0, t1},
     Method -> {"ExplicitRungeKutta", "DifferenceOrder" -> 5},
     PrecisionGoal -> -Log10[relTol], AccuracyGoal -> -Log10[absTol],
     MaxSteps -> Infinity]];

(* ::Section:: *)
(* Check A. Resonance, where the exact answer is known *)
(* y'' + y = cos(t), y(0) = 0, y'(0) = 0 has exact solution y = (t/2) sin(t). *)
fA[t_?NumericQ] := Cos[t];
exact[t_] := t/2 Sin[t];
tols = {10.^-3, 10.^-6, 10.^-9};       (* 1e-3 is MATLAB's default RelTol *)

Print["Check A: resonance, compared with y = (t/2) sin t"];
Print["   RelTol     steps taken    max |error|"];
errCurves = Table[
   Module[{Ys = solveSpring[1., 0., 1., fA, {0., 0.}, {0., 50.}, tol, tol 10.^-3],
           grid, err},
    grid = Flatten[Ys["Grid"]];
    err = Abs[Ys[#][[1]] - exact[#]] & /@ grid;
    Print["   ", sci[tol, 0], "       ", StringPadLeft[ToString[Length[grid] - 1], 8],
      "        ", sci[Max[err]]];
    Transpose[{grid, err + $MachineEpsilon}]],
   {tol, tols}];

plotA = ListLogPlot[errCurves, Joined -> True, PlotRange -> {10^-12, 1},
   Frame -> True, GridLines -> Automatic, FrameLabel -> {"t", "|NDSolve - exact|"},
   PlotLabel -> "Check A: error against the exact solution",
   PlotLegends -> Placed[("RelTol = " <> sci[#, 0] & /@ tols), {0.8, 0.2}],
   ImageSize -> 520];

(* ::Section:: *)
(* Check B. Same problem, two tolerances, long run *)
(* Undamped, unforced: y'' + y = 0, y(0) = 1, y'(0) = 0.
   The exact orbit is the unit circle, forever. *)
fB[t_?NumericQ] := 0.;
YL = solveSpring[1., 0., 1., fB, {1., 0.}, {0., 1000.}, 10.^-3, 10.^-6];   (* MATLAB's defaults *)
YT = solveSpring[1., 0., 1., fB, {1., 0.}, {0., 1000.}, 10.^-8, 10.^-10];  (* tight *)

Print[""];
Print["Check B: undamped spring on [0, 1000], distance from origin at the end"];
Print["   loose tolerances:   ", ToString[NumberForm[Norm[YL[1000.]], {5, 4}]], "   (should be 1)"];
Print["   tight tolerances:   ", ToString[NumberForm[Norm[YT[1000.]], {9, 8}]]];

radius[Ys_] := Transpose[{Flatten[Ys["Grid"]], Norm /@ Ys["ValuesOnGrid"]}];
plotB = ListLinePlot[{radius[YL], radius[YT]}, PlotRange -> {0.95, 1.01},
   Frame -> True, GridLines -> Automatic,
   FrameLabel -> {"t", "distance from origin in the phase plane"},
   PlotLabel -> "Check B: undamped spring, should stay at 1",
   PlotLegends -> Placed[{"RelTol = 1e-3 (MATLAB default)", "RelTol = 1e-8"}, {0.3, 0.2}],
   ImageSize -> 520];

plot = GraphicsRow[{plotA, plotB}, ImageSize -> 1100];
Export["checkSolver.png", plot];
plot                                     (* evaluate on its own line to display *)

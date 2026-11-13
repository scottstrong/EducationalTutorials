(* underTheHood.wl
   Two pictures of what a time stepper is doing. Run as is.
     Part A. Euler and RK4 with the same fixed step, against the exact answer.
     Part B. The step sizes NDSolve actually chooses on a problem with a sudden kick.

   Wolfram port of underTheHood.m.  Numerical only. *)

Quiet@Check[SetDirectory[NotebookDirectory[]], Null];   (* Export writes to Directory[] *)

(* Batch-safe scientific notation, e.g. sci[0.0000521] -> "5.21e-05" *)
sci[x_?NumericQ, d_: 2] := If[x == 0, "0", Module[{e = Floor[Log10[Abs[x]]], mm},
   mm = x/10.^e;
   If[d == 0, ToString[Round[mm]], ToString@NumberForm[mm, {d + 1, d}]] <>
    "e" <> If[e < 0, "-", "+"] <> IntegerString[Abs[e], 10, 2]]];

(* ::Section:: *)
(* Part A. Euler vs. RK4, same step size *)
(* Underdamped spring: y'' + 0.5 y' + y = 0, y(0) = 1, y'(0) = 0. *)
m = 1.;  b = 0.5;  k = 1.;  f[t_] := 0.;
rhs[t_, {yy_, vv_}] := {vv, (f[t] - b vv - k yy)/m};

(* Exact solution of the underdamped free response *)
g  = b/(2 m);                    (* decay rate *)
wd = Sqrt[k/m - g^2];            (* damped frequency *)
y0 = 1.;  v0 = 0.;
exact[t_] := Exp[-g t] (y0 Cos[wd t] + (v0 + g y0)/wd Sin[wd t]);

h = 0.25;                        (* one fixed step size for both methods *)
tt = Range[0., 20., h];
nn = Length[tt];

(* Euler: follow the slope at the start of the step. *)
YE = ConstantArray[0., {nn, 2}];  YE[[1]] = {y0, v0};
Do[YE[[n + 1]] = YE[[n]] + h rhs[tt[[n]], YE[[n]]], {n, 1, nn - 1}];

(* RK4: sample the slope four times across the step and average. *)
YR = ConstantArray[0., {nn, 2}];  YR[[1]] = {y0, v0};
Do[
  k1 = rhs[tt[[n]],       YR[[n]]];
  k2 = rhs[tt[[n]] + h/2, YR[[n]] + h/2 k1];
  k3 = rhs[tt[[n]] + h/2, YR[[n]] + h/2 k2];
  k4 = rhs[tt[[n]] + h,   YR[[n]] + h k3];
  YR[[n + 1]] = YR[[n]] + h/6 (k1 + 2 k2 + 2 k3 + k4),
  {n, 1, nn - 1}];

Print["Part A, step h = ", h, " (", nn - 1, " steps)"];
Print["   Euler max error: ", sci[Max[Abs[YE[[All, 1]] - exact /@ tt]]]];
Print["   RK4   max error: ", sci[Max[Abs[YR[[All, 1]] - exact /@ tt]]]];

plotA = Show[
   Plot[exact[t], {t, 0, 20}, PlotStyle -> {Black, Thick}],
   ListLinePlot[Transpose[{tt, YE[[All, 1]]}], PlotMarkers -> {Automatic, 6}],
   ListPlot[Transpose[{tt, YR[[All, 1]]}], PlotMarkers -> {"\[FilledSquare]", 8},
     PlotStyle -> Red],
   PlotRange -> All, Frame -> True, GridLines -> Automatic,
   FrameLabel -> {"t", "y(t)"}, PlotLabel -> StringForm["Part A: same step h = ``", h],
   Epilog -> Inset[LineLegend[{Black, ColorData[97][1], Red}, {"exact", "Euler", "RK4"}],
     Scaled[{0.85, 0.8}]]];

(* ::Section:: *)
(* Part B. Where NDSolve puts its steps *)
(* Lightly damped spring, released from y = 1 and then pushed at t = 10:
     f(t) = 5 for 10 <= t <= 10.5, and 0 otherwise.
   ?NumericQ hides the formula, so NDSolve cannot see the jumps in advance
   and has to find them by sampling, as ode45 does. *)
mB = 1.;  bB = 0.2;  kB = 1.;
fB[t_?NumericQ] := If[10 <= t <= 10.5, 5., 0.];
rhsB[t_, {yy_, vv_}] := {vv, (fB[t] - bB vv - kB yy)/mB};

YB = Y /. First@NDSolve[
    {Y'[t] == rhsB[t, Y[t]], Y[0] == {1., 0.}},
    Y, {t, 0, 30},
    Method -> {"ExplicitRungeKutta", "DifferenceOrder" -> 5},
    PrecisionGoal -> 6, AccuracyGoal -> 9];
tsteps = Flatten[YB["Grid"]];      (* the times NDSolve actually stepped to *)
hstep  = Differences[tsteps];

Print[""];
Print["Part B: ", Length[hstep], " steps, smallest ",
  sci[Min[hstep]], ", largest ", sci[Max[hstep]]];

plotB1 = ListLinePlot[Transpose[{tsteps, YB["ValuesOnGrid"][[All, 1]]}],
   Mesh -> All, Frame -> True, GridLines -> Automatic, FrameLabel -> {None, "y(t)"},
   PlotLabel -> "Part B: a push at t = 10 (each dot is one step)", ImageSize -> 500,
   AspectRatio -> 0.35];
plotB2 = ListLogPlot[Transpose[{Most[tsteps], hstep}], Joined -> True, Mesh -> All,
   Frame -> True, GridLines -> Automatic, FrameLabel -> {"t", "step size h"},
   ImageSize -> 500, AspectRatio -> 0.35];

plot = GraphicsRow[{plotA, GraphicsColumn[{plotB1, plotB2}]}, ImageSize -> 1100];
Export["underTheHood.png", plot];
plot                                     (* evaluate on its own line to display *)

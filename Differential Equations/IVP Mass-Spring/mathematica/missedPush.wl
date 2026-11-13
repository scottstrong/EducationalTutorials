(* missedPush.wl
   A short push applied to a spring at rest. With its default settings,
   NDSolve can take a step so long that it jumps over the push and never sees it.
   Limiting the step size with MaxStepSize fixes it. Run as is.

   Wolfram port of missedPush.m.  Numerical only.
   The ?NumericQ on f matters: it hides the formula, so NDSolve can only sample
   f, like ode45. Without it NDSolve reads the If, finds the jumps itself, and
   never misses the push.

   What you should see: on this problem NDSolve grows its steps more cautiously
   than ode45 or SciPy and does NOT skip the push (largest |y| about 2.128 vs.
   2.134 with the cap). It got lucky; nothing guarantees it on a shorter or
   later push. MaxStepSize is the guarantee. *)

Quiet@Check[SetDirectory[NotebookDirectory[]], Null];   (* Export writes to Directory[] *)

m = 1.;  b = 0.2;  k = 1.;
f[t_?NumericQ] := If[10 <= t <= 10.5, 5., 0.];   (* push of size 5 for half a second *)
rhs[t_, {yy_, vv_}] := {vv, (f[t] - b vv - k yy)/m};
Y0 = {0., 0.};                                    (* at rest, so nothing happens before t = 10 *)

solve[extra___] := Y /. First@NDSolve[
    {Y'[t] == rhs[t, Y[t]], Y[0] == Y0},
    Y, {t, 0, 30},
    Method -> {"ExplicitRungeKutta", "DifferenceOrder" -> 5},
    PrecisionGoal -> 6, AccuracyGoal -> 9, extra];

YA = solve[];                        (* default step control *)
YB = solve[MaxStepSize -> 0.1];      (* steps capped at 0.1 *)

maxAbsY[Ys_] := Max[Abs[Ys["ValuesOnGrid"][[All, 1]]]];
Print["No MaxStepSize:      largest |y| = ", ToString[NumberForm[maxAbsY[YA], {5, 4}]]];
Print["MaxStepSize = 0.1:   largest |y| = ", ToString[NumberForm[maxAbsY[YB], {5, 4}]]];

curve[Ys_] := Transpose[{Flatten[Ys["Grid"]], Ys["ValuesOnGrid"][[All, 1]]}];
plot = ListLinePlot[{curve[YA], curve[YB]}, PlotRange -> All,
   PlotMarkers -> {{Automatic, 6}, None}, Frame -> True, GridLines -> Automatic,
   FrameLabel -> {"t", "y(t)"}, PlotLabel -> "A push at t = 10, starting from rest",
   PlotLegends -> Placed[{"default step control", "MaxStepSize = 0.1"},
     {Left, Top}], ImageSize -> 1000, AspectRatio -> 0.4];

Export["missedPush.png", plot];
plot                                     (* evaluate on its own line to display *)

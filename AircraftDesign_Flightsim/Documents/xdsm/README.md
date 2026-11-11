# XDSM of the AAMF flight-simulation framework

Draft XDSM diagrams of the computational workflow, generated with
[pyXDSM](https://github.com/mdolab/pyXDSM). Notation follows Lambe & Martins
(2012); style follows the `ex_01_asw` example in `sci-comp-eng-des-opt`.

```bash
pip install pyXDSM
python aamf_xdsm.py             # writes .tikz + .tex (+ .pdf/.png if pdflatex present)
python aamf_xdsm.py --matrix    # text view of the DSM grid, no LaTeX needed
```

| file | contents |
|---|---|
| `aamf_xdsm.py` | the builder — one function per figure, with the modelling rationale in the docstrings |
| `fig_0_aamf_general.{tikz,tex,pdf,png}` | **the framework itself**, no particular aircraft |
| `fig_1_aamf_workflow.{tikz,tex,pdf,png}` | the X-57 chain, published data to CesiumJS |
| `fig_2_force_moment.{tikz,tex,pdf,png}` | inside the force and moment model |
| `fig_3_evaluation.{tikz,tex,pdf,png}` | performance analyses and verification |
| `xdsm_matrix_summary.txt` | text DSM grid, for reviewing before compiling |

**Rendering.** Built with MiKTeX. If `pdflatex` is not on PATH, the script still
writes the TeX sources — `winget install MiKTeX.MiKTeX`, reopen the shell, rerun.
PNG needs no extra tools: PyMuPDF is used when poppler's `pdftoppm` is absent.

---

## Figure 0 — the general case

Read this one first. It answers two questions that the X-57 figures cannot:

**Where does "aircraft-agnostic" actually stop?** Boxes 1–3 are supplied per
aircraft — a data deck, an `Aerodynamics` subclass, and one or more
`PropulsiveElement` subclasses. Boxes 4–16 are framework code that no case
edits. Every example in `examples/` differs only in those three boxes: the F-16
brings `AfterburningTurbojet`, the Cessna brings `PistonEngineModel` +
`PropellerPropulsion`, the 737 brings `CFM56Turbofan`, the X-57 brings a local
`X57ElectricPropulsor`. All five use `CoefficientAerodynamics`.

**One model, five drivers.** The architectural point is the six feedback edges,
five of which land on the same box:

| driver | how it iterates |
|---|---|
| Trim solver | `fmincon` on the residual vector |
| Stability | finite-difference perturbation **sweep**, then `eig` — not a convergence loop |
| Takeoff / landing | time-marching integration, plus a V₁ scan-and-bisect for balanced field length |
| Mission planner | segment loop → node loop → a trim solve per node → fuel integration |
| Simulink S-function | once per solver step inside the 6-DOF time loop |

That is why a defect in the load model surfaces in all five at once.

### Box styles

Shape follows *computational character*, not provenance:

| | |
|---|---|
| `FUNC` | explicit evaluation |
| `IFUNC` | implicit — contains an inner solve |
| `SOLVER` | drives other components to convergence |
| `GROUP` | composite subsystem |

The aircraft-agnostic boundary is marked in the **labels** (“per aircraft”)
rather than by shape, so the two meanings don't collide.

---

## Why separate figures

Figure 0 shows the framework; Figures 1–3 show one aircraft flowing through it.
Merging them would mean losing either the aircraft-agnostic structure in X-57
detail, or the X-57 detail in abstraction.

Within the X-57 set, the force and moment model is evaluated by **two** drivers:
the trim solver, which iterates on it to null five residuals, and the Simulink
S-function, which calls it once per solver step inside the 6-DOF time loop.
Drawing it once and pointing both drivers at it would make the data edges
ambiguous. Figure 1 therefore carries it as a single `GROUP` box with one clean
edge pair per driver, and Figure 2 expands it. Figure 0 does the opposite on
purpose — there the multiplicity of drivers *is* the architecture, so all five
loops are drawn, each labelled with what that driver sends.

## Figure 1 — complete workflow

14 components, 22 couplings, **2 feedback edges**:

* `model → trim` — the trim residual **R**(x, u)
* `sfun → tloop` — forces and mass properties back into the 6-DOF integrator

Everything else is strictly feed-forward, which is what the lower-triangular
structure in the matrix summary shows.

Steps 1–4 build the model (deck → derivatives → configuration → `Aircraft`
object). Steps 5–7 are the trim loop. Steps 8–10 consume the trim point.
Steps 11–13 are the real-time loop. Steps 14–16 are the visualisation path.

## Figure 2 — force and moment model

Given (x, u) this chain is **entirely feed-forward**. There is no internal MDA —
no fixed-point iteration between aerodynamics and propulsion — because
`X57Lookup` recomputes the propeller state itself rather than reading it back
from the propulsive elements. The only iteration inside the model is the
bisection in the propeller disk model (drawn `IFUNC`). That property is what
makes the outer trim solve well-behaved.

**One simplification to be aware of.** The figure shows a single propeller disk
model feeding both the aerodynamic lookup and the propulsive elements. That is
the functional truth, but in the code `X57PropModel` is called *twice* from
different places — once inside `X57Lookup` for the slipstream and thrust
coefficient, once inside each `X57ElectricPropulsor` for the thrust reported to
the load solver. The duplication is deliberate: it makes the aerodynamic
coefficients independent of load-solver evaluation order. The cost is that two
calculations must agree, which is what Phase 3 of the verification suite exists
to enforce.

## Figure 3 — evaluation

No feedback edges: the analyses fan out from the trim point and converge on the
published-value comparison. The iteration is *inside* five of the boxes, drawn
`IFUNC` — blown stall speed is a fixed point (the slipstream depends on the
speed being solved for); cruise shaft power, the level-flight envelope and
V<sub>mc</sub> are each a bisection; and takeoff/landing time-march the ground
roll, rotation, flare and braking phases, with a V₁ scan-and-bisect wrapped
around them for the balanced field length.

---

## Every loop in the framework

Each of these was read out of the classes, not assumed.

| loop | where | mechanism |
|---|---|---|
| **Trim** | `PerformanceAnalysis.solve_trim`, `GenericTrimSolver` | `fmincon`, residual tol 1e-5 |
| **Propeller disk** | `X57PropModel` | bisection, relative 1e-9 on thrust |
| **Stability** | `StabilityAnalysis.linearize` | `for i = 1:nX` perturbation **sweep**, then `eig` — *not* a convergence loop |
| **Takeoff** | `TakeoffAnalysis` | `while` loops for ground roll, rotation, braking; `compute_bfl → find_root_scan → bisect` on V₁ |
| **Landing** | `LandingAnalysis` | `while` loops for approach, flare, rollout |
| **Mission** | `MissionPlanner` | segment loop → node loop → trim per node → fuel integration; no outer weight loop |
| **6-DOF time** | Simulink variable-step solver | per-step, not to a residual |
| **Performance inner solves** | `X57Performance` | fixed point (blown stall); bisection (cruise power, V<sub>H</sub>, V<sub>mc</sub>) |

The X-57 case exercises all of these except `MissionPlanner`, which `X57Simp`
does not call — which is why it appears in Figure 0 but not Figure 1.

### On the trim being five residuals, not six

Cruise solves (α, δ<sub>e</sub>, throttle, δ<sub>a</sub>, δ<sub>r</sub>) at
γ = 0. The DiTTo case fixes throttle at the published 100 % torque setting and
solves for γ instead, because the aircraft climbs at that speed and power rather
than flying level.

Residuals are (u̇, ẇ, ṗ, q̇, ṙ) — side force is deliberately excluded. The
single left seat puts the CG 23 mm off centreline; that can be trimmed in roll
with a little aileron, but the resulting side force has nothing to balance
against in steady wings-level flight. It settles at about 8 N, 0.06 % of weight.
Claiming six converged residuals would mean claiming a solution that does not
exist.

### On the real-time loop being open in controls

The 6-DOF loop closes through the **state** only. Controls come open-loop from a
From Workspace block holding **u**\*. The net_ctrl return path from CesiumJS
exists in the shared Simulink model — UDP Receive and Unpack blocks are both
present — but the Unpack block has **zero output ports enabled**, so nothing
flows back. The diagram shows the forward path only, which is what the model
actually does.

---

## Known gaps this diagram does not paper over

* **φ/θ transposition.** Real net_fdm packets carry pitch in the `phi` field and
  roll in `theta`. This is in the shared `FlightSim_Aircraft.slx`, so it affects
  every example, and is not yet fixed. Figure 1's `fdm` box is therefore correct
  as *drawn* but the model behind it transposes two fields.
* **Packet rate.** The net_fdm packer inherits the solver step rather than
  running on a fixed sample time, producing ~865 packets/s against a viewer that
  reads at 10 Hz.

## Sources

* Lambe, A. B. and Martins, J. R. R. A., "Extensions to the Design Structure
  Matrix for the Description of Multidisciplinary Design, Analysis, and
  Optimization Processes", *Structural and Multidisciplinary Optimization*,
  46(2), 2012.
* Martins, J. R. R. A. and Ning, A., *Engineering Design Optimization*,
  Cambridge University Press, 2021.
* Model sources: see the header of
  `AircraftDesign_Flightsim/MATLAB/examples/X57/X57Data.m`.

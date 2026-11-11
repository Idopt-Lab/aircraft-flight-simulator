"""pyXDSM diagrams of the AAMF flight-simulation framework.

Renders four linked XDSM figures - one general, three for the X-57 case:

    fig_0_aamf_general      the framework itself, no particular aircraft
    fig_1_aamf_workflow     the X-57 chain, data deck to CesiumJS
    fig_2_force_moment      inside the force and moment model
    fig_3_evaluation        performance analyses and the verification suite

Figure 0 is the one to read first: it shows where the aircraft-agnostic line
falls (boxes 1-3 are supplied per aircraft, 4-16 are framework code) and that
there is ONE force and moment model with FIVE drivers iterating on it.
Figures 1-3 are the X-57 instance of that architecture.

Run from anywhere::

    python aamf_xdsm.py

pyXDSM writes .tikz and .tex unconditionally; producing the PDF needs
``pdflatex`` (TeX Live or MiKTeX) on PATH, and the PNG additionally needs
``pdftoppm`` (poppler).  If those are missing the script still writes the
TeX sources and says so, so the diagrams can be rendered elsewhere.

--------------------------------------------------------------------------
WHY SEPARATE FIGURES AND NOT ONE
--------------------------------------------------------------------------
Figure 0 shows the framework; Figures 1-3 show one aircraft flowing through
it.  Merging them would mean either losing the aircraft-agnostic structure
in X-57 detail, or losing the X-57 detail in abstraction.

Within the X-57 set, the force and moment model is evaluated by TWO
different drivers: the trim solver, which iterates on it to null five
residuals, and the Simulink S-function, which calls it once per solver step
inside the 6-DOF time loop.  Drawing the model once and letting both drivers
point at it makes the data edges ambiguous - a given edge would mean
different things depending on which loop is active.  So Figure 1 carries the
model as a single GROUP box with one clean edge pair per driver, and Figure
2 expands it.

Figure 0 does the opposite on purpose: it draws ALL five drivers with their
loops back into the one shared model, because there the multiplicity IS the
architecture, and each edge is labelled with what that driver sends.

STYLE.  Box shape follows computational character, not provenance:
  FUNC    explicit evaluation            IFUNC   implicit, has an inner solve
  SOLVER  drives others to convergence   GROUP   composite subsystem
The aircraft-agnostic boundary is marked in the LABELS instead - boxes 1-3
of Figure 0 read "per aircraft"; everything else is framework code.

--------------------------------------------------------------------------
THE ITERATIVE STRUCTURE, WHICH IS THE POINT OF THE DIAGRAMS
--------------------------------------------------------------------------
Every one of these was read out of the classes rather than assumed.

In the FRAMEWORK generally (Fig 0), five drivers iterate on the one shared
force and moment model:

  A  TRIM                PerformanceAnalysis.solve_trim / GenericTrimSolver,
                         fmincon against a residual tolerance.
  B  STABILITY           StabilityAnalysis.linearize runs a finite-difference
                         perturbation SWEEP (for i = 1:nX) and then eig.  A
                         sweep, not a convergence loop - which is why it is
                         drawn FUNC and not SOLVER.
  C  TAKEOFF / LANDING   Time-marching integration: while-loops for ground
                         roll, rotation and braking in TakeoffAnalysis, and
                         for approach, flare and rollout in LandingAnalysis.
                         On a multi-engine aircraft calculate_takeoff
                         dispatches to compute_bfl, which scans and then
                         bisects on V1 for the balanced field length - a root
                         find wrapped AROUND the integrations.
  D  MISSION             MissionPlanner loops over segments, then over nodes,
                         solving a trim at each node, then integrates fuel.
                         No outer weight-convergence loop.
  E  SIMULINK            The 6-DOF block integrates what the S-function
                         returns and feeds the state back next step.

In the X-57 CASE (Figs 1-3), four more specific loops:

  1  TRIM (Fig 1, Fig 2).  Five variables against five residuals.  Cruise
     solves (alpha, delta_e, throttle, delta_a, delta_r) at gamma = 0; the
     DiTTo case fixes throttle at the published 100% torque setting and
     solves for gamma instead, because the aircraft climbs at that speed and
     power rather than flying level.  Residuals are (udot, wdot, pdot, qdot,
     rdot) - five, not six.  Side force is deliberately excluded: the single
     left seat puts the CG 23 mm off centreline, which can be trimmed in roll
     but leaves a residual Fy that steady wings-level flight has nothing to
     balance against.

  2  PROPELLER DISK (Fig 2).  X57PropModel inverts actuator-disk momentum
     theory by bisection for the thrust a commanded shaft power produces.
     Drawn IFUNC because it is implicit.

  3  6-DOF TIME LOOP (Fig 1).  Controls do NOT close this loop: they come
     open-loop from a From Workspace block holding u*.  The net_ctrl return
     path from CesiumJS exists in the model but its Unpack block has zero
     output ports enabled, so nothing flows back.

  4  PERFORMANCE INNER SOLVES (Fig 3).  Blown stall speed is a fixed point
     (the slipstream depends on the speed being solved for); cruise shaft
     power, level-flight envelope and Vmc are each a bisection; takeoff and
     landing are the time-marching integrations of C above.

Worth stating explicitly because it is what makes the trim tractable: given
(x, u), the force and moment model in Figure 2 is entirely FEED-FORWARD.
There is no internal MDA - no fixed-point iteration between aerodynamics
and propulsion - because X57Lookup recomputes the propeller state itself
rather than reading it back from the propulsive elements.  The only
iteration inside the model is the disk-model bisection.

--------------------------------------------------------------------------
A MODELLING CHOICE THE DIAGRAM MAKES LOOK CLEANER THAN THE CODE
--------------------------------------------------------------------------
Figure 2 shows ONE propeller disk model feeding both the aerodynamic lookup
and the propulsive elements.  That is the functional truth, but in the code
X57PropModel is called TWICE from different places - once inside X57Lookup
for the slipstream and thrust coefficient, once inside each
X57ElectricPropulsor for the thrust reported to the load solver.  The
duplication is deliberate: it makes the aerodynamic coefficients
independent of the order the framework evaluates its load solvers in.  The
cost is that two calculations must agree, and Phase 3 of the verification
suite exists mainly to enforce that.

--------------------------------------------------------------------------
SOURCES FOR THE NOTATION
--------------------------------------------------------------------------
  Lambe, A.B. and Martins, J.R.R.A., "Extensions to the Design Structure
    Matrix for the Description of Multidisciplinary Design, Analysis, and
    Optimization Processes", Structural and Multidisciplinary
    Optimization, 46(2), 2012 - the XDSM convention itself.
  Martins, J.R.R.A. and Ning, A., "Engineering Design Optimization",
    Cambridge University Press, 2021 - process/data line conventions.
  Style follows sci-comp-eng-des-opt/src/aircraft_sizing/examples/ex_01_asw
    /viz/xdsm/asw_sizing_xdsm.py.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

OUTPUT_DIR = Path(__file__).resolve().parent
SCRIPT_COMMAND = f"python {Path(__file__).name}"


class MissingPrerequisiteError(RuntimeError):
    """Raised when rendering prerequisites are unavailable."""


def _pyxdsm_api():
    try:
        from pyxdsm.XDSM import FUNC, GROUP, IFUNC, RIGHT, SOLVER, XDSM
    except ModuleNotFoundError as exc:  # pragma: no cover
        raise MissingPrerequisiteError(
            "Missing prerequisite: Python package 'pyXDSM' is not installed in the "
            f"active environment. Run `pip install pyXDSM`, then rerun `{SCRIPT_COMMAND}`."
        ) from exc
    return XDSM, SOLVER, FUNC, IFUNC, GROUP, RIGHT


# =========================================================================
#  FIGURE 0 - the general AAMF architecture, no particular aircraft
# =========================================================================
def build_general_xdsm():
    """AAMF in general: what the framework provides and what a case supplies.

    The framework is aircraft-agnostic, and this figure is about where that
    line falls.  A new aircraft supplies four things - a data deck, an
    Aerodynamics implementation, one or more PropulsiveElement
    implementations, and a build script.  Everything downstream of the
    assembly is generic and identical for every case in examples/.  The two
    pluggable boxes are drawn GROUP to mark that boundary; every FUNC and
    SOLVER box is framework code that a case never edits.

    The architectural point the diagram is meant to make is that there is
    ONE force and moment model and FIVE drivers that iterate on it:

        trim solver         fmincon on 5 or 6 residuals
        stability           finite-difference perturbation sweep
        takeoff / landing   time-marching integration, plus a V1 root find
                            for balanced field length on multi-engine cases
        mission planner     segment loop, node loop, a trim solve per node
        Simulink S-function once per solver step inside the 6-DOF time loop

    That is why every one of them has an edge back into the load model, and
    why a defect in the model shows up in all five at once.

    Verified against the classes rather than assumed:
      PerformanceAnalysis / GenericTrimSolver  fmincon, residual tolerance
      StabilityAnalysis.linearize              for i = 1:nX perturbation
                                               sweep, then eig - a sweep,
                                               NOT a convergence loop
      TakeoffAnalysis                          while-loops for ground roll,
                                               rotation and braking;
                                               compute_bfl -> find_root_scan
                                               -> bisect on V1
      LandingAnalysis                          while-loops for approach,
                                               flare and rollout
      MissionPlanner                           for over segments, for over
                                               nodes, trim per node, then
                                               fuel integration
    """
    XDSM, SOLVER, FUNC, IFUNC, GROUP, RIGHT = _pyxdsm_api()
    x = XDSM(use_sfmath=False)

    # ---- case-supplied layer (GROUP marks the plug-in boundary) --------
    x.add_system("deck", FUNC, (r"\text{1: Aircraft data deck}", r"\textit{per aircraft}"))
    x.add_system(
        "aeromod",
        GROUP,
        (r"\text{2: Aerodynamics model}", r"\texttt{Aerodynamics}\ \text{subclass}",
         r"\textit{per aircraft}"),
    )
    x.add_system(
        "propmod",
        GROUP,
        (r"\text{3: Propulsion models}", r"\texttt{PropulsiveElement}\ \text{subclass}",
         r"\textit{per aircraft}"),
        stack=True,
    )

    # ---- framework core ------------------------------------------------
    x.add_system("cfg", FUNC, (r"\text{4: Aircraft assembly}", r"\texttt{AircraftConfigurator}"))
    x.add_system("trim", SOLVER, (r"\text{5, 7: Trim solver}", r"\texttt{fmincon} \text{ on } \mathcal{R}"))
    x.add_system(
        "loads",
        GROUP,
        (r"\text{6: Force and moment}", r"\text{model, shared}"),
    )

    # ---- analyses, each re-entering the load model ---------------------
    x.add_system("perfpt", FUNC, (r"\text{8: Performance point}", r"\texttt{evaluate\_trim\_point}"))
    x.add_system("stab", FUNC, (r"\text{9: Stability}", r"\text{FD sweep} \rightarrow \mathbf{A},\mathbf{B},\ \text{eig}"))
    x.add_system("field", IFUNC, (r"\text{10: Takeoff / landing}", r"\text{time march} + V_1"))
    x.add_system("mission", SOLVER, (r"\text{11: Mission planner}", r"\text{segment, node, trim}"))

    # ---- export and real-time path -------------------------------------
    x.add_system("xport", FUNC, (r"\text{12: Simulink export}", r"\text{base workspace}"))
    x.add_system("sim6", SOLVER, (r"\text{13, 15: 6-DOF}", r"\text{time integration}"))
    x.add_system("sfun", FUNC, (r"\text{14: S-function}", r"\texttt{aircraft\_sfunc}"))
    x.add_system("viz", FUNC, (r"\text{16: net\_fdm, UDP,}", r"\text{CesiumJS}"))

    # ---- external inputs ------------------------------------------------
    x.add_input("deck", (r"\text{geometry, mass, limits,}", r"\text{published or test data}"))
    x.add_input("aeromod", (r"\text{coefficient tables,}", r"\text{DATCOM, or callback}"))
    x.add_input("propmod", (r"\text{engine deck, motor}", r"\text{or propeller ratings}"))
    x.add_input("trim", (r"\text{condition, variables,}", r"\text{bounds, tolerance}"))
    x.add_input("loads", (r"\text{ISA atmosphere}", r"\texttt{isa1976}"))
    x.add_input("field", (r"\text{runway}, h, \mu,", r"\Delta t_{field}\ \text{(own)}"))
    x.add_input("mission", r"\text{segment definitions}")
    x.add_input("sim6", (r"t_{end},\ \Delta t_{sim}", r"\text{(Simulink only)}"))

    # ---- definition, feed-forward ---------------------------------------
    x.connect("deck", "aeromod", r"S,b,\bar c,\ \text{polars}")
    x.connect("deck", "propmod", r"\text{ratings, installation}")
    x.connect("deck", "cfg", (r"\text{components, stations,}", r"\text{surface limits}"))
    x.connect("aeromod", "cfg", r"\texttt{Aerodynamics}")
    x.connect("propmod", "cfg", r"\texttt{PropulsiveElement}\ \text{objects}")
    x.connect("cfg", "trim", r"m,\ \mathbf{r}_{cg},\ \mathbf{I}")
    x.connect("cfg", "loads", (r"\texttt{ac}\text{: frames,}", r"\text{components, solvers}"))
    x.connect("cfg", "xport", r"\texttt{ac},\ n_{cs},n_{pe}")

    # ---- the five drivers of the one shared model ------------------------
    x.connect("trim", "loads", r"\mathbf{x}^{(k)},\ \mathbf{u}^{(k)}")
    x.connect("loads", "trim", r"\mathcal{R}")                      # FEEDBACK

    x.connect("trim", "perfpt", r"\mathbf{x}^*,\ \mathbf{u}^*")
    x.connect("trim", "stab", r"\mathbf{x}^*,\ \mathbf{u}^*")
    x.connect("trim", "field", r"\mathbf{x}^*,\ \mathbf{u}^*")
    x.connect("trim", "mission", r"\mathbf{x}^*,\ \mathbf{u}^*")
    x.connect("trim", "xport", r"\mathbf{x}^*,\ \mathbf{u}^*")

    x.connect("loads", "perfpt", r"\mathbf{F},\mathbf{M},\ \dot m_{fuel}")
    x.connect("perfpt", "loads", (r"\mathbf{x}^*,\mathbf{u}^*\ \text{once}", r"\text{(no iteration)}"))
    x.connect("loads", "stab", r"\mathbf{F},\mathbf{M}")
    x.connect("stab", "loads", (r"\mathbf{x}\pm\delta,\ \mathbf{u}\pm\delta", r"\text{FD sweep}"))   # FEEDBACK
    x.connect("loads", "field", (r"T\ \text{(Propulsion solver)},", r"C_L,C_D\ \text{(Aero solver)}"))
    x.connect("field", "loads", (r"\mathbf{x}(t),\ \Delta t = 0.05\ \text{s}", r"\text{own ground EOM}"))  # FEEDBACK
    x.connect("mission", "trim", (r"\text{node condition}", r"\text{per segment}"))  # FEEDBACK

    # ---- export and the real-time loop -----------------------------------
    x.connect("xport", "sim6", r"\text{initial state}")
    x.connect("xport", "sfun", r"\texttt{ac},\ n_{total},\ \mathbf{u}(t)")
    x.connect("sim6", "sfun", r"\mathbf{x}(t)")
    x.connect("sfun", "sim6", r"\mathbf{F},\mathbf{M},m,\mathbf{I}")     # FEEDBACK
    x.connect("loads", "sfun", r"\mathbf{F},\mathbf{M}")
    x.connect("sfun", "loads", r"\mathbf{x},\mathbf{u}\ \text{per step}")  # FEEDBACK
    x.connect("sim6", "viz", r"\mathbf{r}_{NED},\ \phi,\theta,\psi")

    # ---- deliverables -----------------------------------------------------
    x.add_output("trim", r"\mathbf{x}^*,\ \mathbf{u}^*", side=RIGHT)
    x.add_output("perfpt", r"L/D,\ P,\ \text{margins}", side=RIGHT)
    x.add_output("stab", r"\mathbf{A},\mathbf{B},\ \lambda_i,\ \text{modes}", side=RIGHT)
    x.add_output("field", r"s_{TO},\ s_{LDG},\ V_1", side=RIGHT)
    x.add_output("mission", r"\text{fuel, time, range}", side=RIGHT)
    x.add_output("viz", r"\text{live 3-D view}", side=RIGHT)

    # ---- process: the main sequence plus every loop ------------------------
    x.add_process(["deck", "aeromod", "cfg", "trim"], arrow=True)
    x.add_process(["deck", "propmod", "cfg"], arrow=True)
    x.add_process(["trim", "loads", "trim"], arrow=True)              # trim loop
    x.add_process(["trim", "perfpt", "stab"], arrow=True)
    x.add_process(["stab", "loads", "stab"], arrow=True)              # FD sweep
    x.add_process(["stab", "field"], arrow=True)
    x.add_process(["field", "loads", "field"], arrow=True)            # time march
    x.add_process(["field", "mission"], arrow=True)
    x.add_process(["mission", "trim", "mission"], arrow=True)         # node loop
    x.add_process(["mission", "xport", "sim6"], arrow=True)
    x.add_process(["sim6", "sfun", "sim6"], arrow=True)               # time loop
    x.add_process(["sfun", "loads", "sfun"], arrow=True)              # per step
    x.add_process(["sim6", "viz"], arrow=True)

    return x


# =========================================================================
#  FIGURE 1 - the complete framework
# =========================================================================
def build_workflow_xdsm():
    """Complete chain: published data through to the CesiumJS viewer.

    Diagonal order is chosen so that the only entries below the diagonal are
    the two genuine feedback edges - the trim residual returning to the trim
    solver, and the forces returning to the 6-DOF integrator.
    """
    XDSM, SOLVER, FUNC, IFUNC, GROUP, RIGHT = _pyxdsm_api()
    x = XDSM(use_sfmath=False)

    # ---- model definition, feed-forward -------------------------------
    x.add_system("deck", FUNC, (r"\text{1: Input deck}", r"\texttt{X57Data}"))
    x.add_system("derv", FUNC, (r"\text{2: Derivative build-up}", r"\texttt{X57Aero}"))
    x.add_system("ctx", FUNC, (r"\text{3: Configuration}", r"\texttt{X57Context}"))
    x.add_system("asm", FUNC, (r"\text{4: Aircraft assembly}", r"\texttt{X57Build}"))

    # ---- trim, iterative ----------------------------------------------
    x.add_system("trim", SOLVER, (r"\text{5, 7: Trim solver}", r"\text{5 vars, 5 residuals}"))
    x.add_system(
        "model",
        GROUP,
        (r"\text{6: Force and moment}", r"\text{model (Fig. 2)}"),
    )

    # ---- post-trim analyses -------------------------------------------
    x.add_system("stab", FUNC, (r"\text{8: Linearisation}", r"\text{FD sweep} \rightarrow \mathbf{A},\mathbf{B}"))
    x.add_system("perf", FUNC, (r"\text{9: Performance and}", r"\text{field (Fig. 3)}"))

    # ---- export and real-time simulation ------------------------------
    x.add_system("xport", FUNC, (r"\text{10: Simulink export}", r"\text{base workspace}"))
    x.add_system("tloop", SOLVER, (r"\text{11, 13: 6-DOF}", r"\text{time integration}"))
    x.add_system("sfun", FUNC, (r"\text{12: S-function}", r"\texttt{aircraft\_sfunc}"))
    x.add_system("fdm", FUNC, (r"\text{14: Geodetic and}", r"\text{net\_fdm packing}"))
    x.add_system("bridge", FUNC, (r"\text{15: UDP/WebSocket}", r"\texttt{main.py}"))
    x.add_system("view", FUNC, (r"\text{16: CesiumJS viewer}", r"\texttt{main.js}"))

    # ---- external inputs ----------------------------------------------
    x.add_input("deck", (r"\text{published X-57 data:}", r"\text{NTRS, ICAS, ASTM F44}"))
    x.add_input("ctx", (r"\text{gear},\ N_{CM},\ N_{HLP},", r"\text{windmilling}"))
    x.add_input("trim", (r"V,\ h,\ \text{case}", r"\mathbf{u}^{(0)},\ \mathbf{u}_{lb},\mathbf{u}_{ub}"))
    x.add_input("tloop", r"t_{end},\ \Delta t")
    x.add_input("bridge", r"\text{ports } 50000/8765")

    # ---- feed-forward couplings ---------------------------------------
    x.connect("deck", "derv", (r"S,b,\bar c,\lambda,AR,", r"\text{airfoil}, S_h, S_v, l_h"))
    x.connect("deck", "ctx", r"\text{deck}")
    x.connect("deck", "asm", (r"\text{mass table, stations,}", r"\text{control limits}"))
    x.connect("deck", "perf", r"\text{published targets}")
    x.connect("derv", "ctx", (r"C_{L\alpha},C_{m\alpha},C_{m0},\varepsilon_0,", r"\text{polars, blown lift}"))
    x.connect("ctx", "asm", r"\texttt{ctx}")
    x.connect("ctx", "model", (r"\mathbf{a},\ N_{CM},\ N_{HLP},", r"\text{gear}"))
    x.connect("asm", "trim", r"m,\ \mathbf{r}_{cg},\ \mathbf{I}")
    x.connect("asm", "model", r"\texttt{ac}\ \text{(frames, elements)}")
    x.connect("asm", "xport", r"\texttt{ac},\ n_{cs},n_{pe}")

    x.connect("trim", "model", r"\mathbf{x}^{(k)},\ \mathbf{u}^{(k)}")
    x.connect("model", "trim", r"\mathcal{R}(\mathbf{x},\mathbf{u})")  # FEEDBACK

    x.connect("trim", "stab", r"\mathbf{x}^*,\ \mathbf{u}^*")
    x.connect("trim", "perf", r"\mathbf{x}^*,\ \mathbf{u}^*")
    x.connect("trim", "xport", r"\mathbf{x}^*,\ \mathbf{u}^*")

    x.connect("xport", "tloop", (r"\text{Initialpos, Vel,}", r"\text{Ori, Rot}"))
    x.connect("xport", "sfun", (r"\texttt{ac},\ n_{total},", r"\texttt{control\_input\_data} \rightarrow \mathbf{u}"))

    x.connect("tloop", "sfun", r"\mathbf{x}(t)\ (12)")
    x.connect("sfun", "tloop", r"\mathbf{F},\mathbf{M},m,\mathbf{I}")  # FEEDBACK
    x.connect("tloop", "fdm", r"\mathbf{r}_{NED},\ \phi,\theta,\psi")
    x.connect("fdm", "bridge", (r"\text{net\_fdm v24,}", r"\text{408 B, UDP}"))
    x.connect("bridge", "view", r"\text{JSON / WebSocket}")

    # ---- deliverables --------------------------------------------------
    x.add_output("derv", r"x_{np},\ \text{SM},\ \mathbf{a}", side=RIGHT)
    x.add_output("trim", r"\mathbf{x}^*,\ \mathbf{u}^*", side=RIGHT)
    x.add_output("stab", r"\mathbf{A},\mathbf{B},\ \lambda_i", side=RIGHT)
    x.add_output("perf", r"V_s,\ L/D,\ V_{mc},\ RC,\ V_H", side=RIGHT)
    x.add_output("view", r"\text{rendered trajectory}", side=RIGHT)

    # ---- process flow ---------------------------------------------------
    x.add_process(["deck", "derv", "ctx", "asm", "trim"], arrow=True)
    x.add_process(["trim", "model", "trim"], arrow=True)          # trim loop
    x.add_process(["trim", "stab", "perf", "xport", "tloop"], arrow=True)
    x.add_process(["tloop", "sfun", "tloop"], arrow=True)          # time loop
    x.add_process(["tloop", "fdm", "bridge", "view"], arrow=True)

    return x


# =========================================================================
#  FIGURE 2 - inside the force and moment model
# =========================================================================
def build_force_moment_xdsm():
    """One evaluation of the force and moment model at a given state.

    Called by the trim solver during trim and by the S-function at every
    Simulink step.  Given (x, u) the whole chain is feed-forward; the only
    iteration is the bisection inside the propeller disk model.
    """
    XDSM, SOLVER, FUNC, IFUNC, GROUP, RIGHT = _pyxdsm_api()
    x = XDSM(use_sfmath=False)

    x.add_system("call", SOLVER, (r"\text{1, 6: Trim solver}", r"\text{or S-function}"))
    x.add_system("cond", FUNC, (r"\text{2: Flight condition}", r"\text{and atmosphere}"))
    x.add_system(
        "disk",
        IFUNC,
        (r"\text{3: Propeller disk}", r"\texttt{X57PropModel}", r"\text{bisection}"),
        stack=True,
    )
    x.add_system("aero", FUNC, (r"\text{4: Aerodynamic}", r"\text{coefficients}", r"\texttt{X57Lookup}"))
    x.add_system(
        "elec",
        FUNC,
        (r"\text{4: Propulsive elements}", r"\texttt{X57ElectricPropulsor}"),
        stack=True,
    )
    x.add_system("grav", FUNC, (r"\text{4: Gravity}", r"\texttt{GravityLoadSolver}"))
    x.add_system("sum", FUNC, (r"\text{5: Load assembly}", r"\text{and residual}"))

    x.add_input("cond", r"h,\ \text{ISA}")
    # Note: a label line must not START with "[" - inside an array, "\\ [" is
    # parsed as the optional vertical-space argument of \\ and LaTeX then tries
    # to read the contents as a length.
    x.add_input("disk", (r"D,\ \eta_{FoM},", r"J_{des},\ J_{0}\ \text{(fixed pitch)}"))
    x.add_input("aero", (r"\mathbf{a}\ \text{(Fig. 1)},", r"\text{gear, flap}"))
    x.add_input("elec", r"Q_{max},\ N_{cmd},\ P_{max}")
    x.add_input("grav", r"m,\ \mathbf{r}_{cg}")

    x.connect("call", "cond", r"\mathbf{x},\ \mathbf{u}")
    x.connect("cond", "disk", (r"V_{ax}=u,\ \rho", r"\text{(body-x, not } |V| )"))
    x.connect("cond", "aero", (r"\alpha,\beta,\bar q,M,", r"\hat p,\hat q,\hat r"))
    x.connect("cond", "elec", r"\rho,\ V_{ax}")
    x.connect("cond", "grav", r"\phi,\theta,\psi")

    x.connect("disk", "aero", (r"(V_{slip}/V)^2 \rightarrow \text{blow}", r"T_{CM} \rightarrow T_c \rightarrow K_0,K_2"))
    x.connect("disk", "elec", r"T_{each},\ P_{ideal}")

    x.connect("aero", "sum", r"C_L,C_D,C_Y,C_l,C_m,C_n")
    x.connect("elec", "sum", (r"\mathbf{F}_T,\ \mathbf{M}_T,", r"P_{shaft},P_{elec}"))
    x.connect("grav", "sum", r"\mathbf{W}_b")

    x.connect("sum", "call", (r"\mathcal{R}=(\dot u,\dot w,\dot p,\dot q,\dot r)", r"\text{5 of 6}"))  # FEEDBACK

    x.add_output("sum", (r"\mathbf{F},\ \mathbf{M}\ \text{about cg},", r"\dot m_{fuel}=0"), side=RIGHT)
    x.add_output("elec", r"P_{elec},\ I_{bus}", side=RIGHT)
    x.add_output("aero", r"\text{envelope flags}", side=RIGHT)

    x.add_process(["call", "cond", "disk", "aero", "sum", "call"], arrow=True)
    x.add_process(["disk", "elec", "sum"], arrow=True)
    x.add_process(["cond", "grav", "sum"], arrow=True)

    return x


# =========================================================================
#  FIGURE 3 - evaluation: performance analyses and verification
# =========================================================================
def build_evaluation_xdsm():
    """What is done with the model once it trims.

    Four of these analyses carry their own inner solve, drawn IFUNC.  The
    verification suite is included because it is part of the research
    workflow, not an afterthought: it consumes the deck, the model, the trim
    point, the Simulink base workspace and the decoded UDP packets.
    """
    XDSM, SOLVER, FUNC, IFUNC, GROUP, RIGHT = _pyxdsm_api()
    x = XDSM(use_sfmath=False)

    x.add_system("pt", FUNC, (r"\text{1: Trim point and}", r"\text{configuration}"))
    x.add_system(
        "stall",
        IFUNC,
        (r"\text{2: Minimum steady}", r"\text{flight speed}", r"\text{fixed point}"),
    )
    x.add_system("ld", FUNC, (r"\text{3: Lift-to-drag}", r"\text{sweep}"))
    x.add_system("pwr", IFUNC, (r"\text{4: Cruise shaft power}", r"\text{disk model inverted}"))
    x.add_system("env", IFUNC, (r"\text{5: Level-flight}", r"\text{envelope } V_H(h)"))
    x.add_system("vmc", IFUNC, (r"\text{6: } V_{mc}\ \text{one motor out}", r"\text{rudder-limited}"))
    x.add_system("rc", FUNC, (r"\text{7: DiTTo climb}", r"RC=(T-D)V/W"))
    # TakeoffAnalysis and LandingAnalysis are NOT plain evaluations.  Both
    # time-march the ground roll, rotation and braking phases, and because the
    # X-57 has more than one propulsive element calculate_takeoff dispatches to
    # compute_bfl, which scans and then bisects on V1 for the balanced field
    # length - a root find wrapped around the integrations.
    x.add_system(
        "fld",
        IFUNC,
        (r"\text{8: Takeoff and landing}", r"\text{time march} + V_1\ \text{bisection}"),
    )
    x.add_system("cmp", FUNC, (r"\text{9: Comparison vs}", r"\text{published targets}"))
    x.add_system(
        "ver",
        FUNC,
        (r"\text{10: Verification suite}", r"\text{6 phases}"),
        stack=True,
    )

    x.add_input("pt", (r"\mathbf{x}^*,\mathbf{u}^*\ \text{(Fig. 1)},", r"\texttt{ac},\ \mathbf{a}"))
    x.add_input("stall", r"C_{L max}(\delta_f)")
    x.add_input("pwr", r"\eta_{motor},\ E_{batt}")
    x.add_input("vmc", r"\delta_{r,max},\ y_{eng}")
    x.add_input("cmp", r"\texttt{d.targets}")
    x.add_input("ver", (r"\text{deck, model, trim,}", r"\text{.slx, net\_fdm}"))

    x.connect("pt", "stall", r"W,\ \rho_0,\ S")
    x.connect("pt", "ld", r"\text{polars},\ h_{cr}")
    x.connect("pt", "pwr", r"D_{cr},\ V_{cr}")
    x.connect("pt", "env", r"P_{max},\ \text{polar}")
    x.connect("pt", "vmc", r"C_{n\delta r},\ T(V)")
    x.connect("pt", "rc", r"T_{total},\ D")
    x.connect("pt", "fld", r"\texttt{takeoff\_params}")

    x.connect("stall", "rc", r"V_{s,blown} \rightarrow V_y")

    x.connect("stall", "cmp", r"V_s,V_{s1},V_{s0},V_{s,bl}")
    x.connect("ld", "cmp", r"(L/D)_{max},\ (L/D)_{cr}")
    x.connect("pwr", "cmp", r"P_{shaft},\ \eta_p,\ E,\ R")
    x.connect("env", "cmp", r"V_H(h)")
    x.connect("vmc", "cmp", r"V_{mc}")
    x.connect("rc", "cmp", r"RC,\ I_{bus}")
    x.connect("fld", "cmp", r"s_{TO},\ s_{LDG}")
    x.connect("cmp", "ver", r"\text{model vs published}")

    x.add_output("cmp", (r"\text{comparison table,}", r"\texttt{X57Performance.csv}"), side=RIGHT)
    x.add_output("ver", (r"\text{GREEN/RED per check,}", r"\texttt{all\_pass}"), side=RIGHT)

    x.add_process(
        ["pt", "stall", "ld", "pwr", "env", "vmc", "rc", "fld", "cmp", "ver"], arrow=True
    )

    return x


# =========================================================================
#  Rendering
# =========================================================================
DIAGRAMS = (
    ("fig_0_aamf_general", build_general_xdsm, "General AAMF architecture"),
    ("fig_1_aamf_workflow", build_workflow_xdsm, "Complete AAMF workflow"),
    ("fig_2_force_moment", build_force_moment_xdsm, "Force and moment model"),
    ("fig_3_evaluation", build_evaluation_xdsm, "Performance and verification"),
)


def _cleanup_latex_artifacts(stem: Path) -> None:
    for suffix in (".aux", ".fdb_latexmk", ".fls", ".log"):
        artifact = stem.with_suffix(suffix)
        if artifact.exists():
            artifact.unlink()


def _render_pdf(tex_path: Path) -> Path | None:
    """Compile one figure.  A failure is reported, not raised.

    One figure with a bad label should not stop the others from building, and
    the useful part of a LaTeX failure is the few lines around the error, not
    the Python traceback - so those are echoed from the .log instead.
    """
    if shutil.which("pdflatex") is None:
        return None

    result = subprocess.run(
        ["pdflatex", "-halt-on-error", "-interaction=nonstopmode", tex_path.name],
        cwd=tex_path.parent,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print(f"    pdflatex FAILED for {tex_path.name}:", file=sys.stderr)
        for line in result.stdout.splitlines():
            if line.startswith("!") or line.startswith("l."):
                print(f"      {line}", file=sys.stderr)
        return None

    _cleanup_latex_artifacts(tex_path.with_suffix(""))
    return tex_path.with_suffix(".pdf")


PNG_DPI = 300


def _render_png(pdf_path: Path) -> Path | None:
    """Rasterise the PDF, preferring poppler and falling back to PyMuPDF.

    PyMuPDF is a self-contained wheel, so on Windows it avoids having to
    install poppler or Ghostscript just to get a picture.  PNG rather than
    JPEG on purpose: these are line drawings with text, and JPEG's chroma
    subsampling puts visible ringing around every glyph and rule.  If a JPEG
    is genuinely needed for a template that refuses PNG, convert from the
    PNG afterwards rather than rasterising to JPEG directly.
    """
    png_path = pdf_path.with_suffix(".png")

    pdftoppm = shutil.which("pdftoppm")
    if pdftoppm is not None:
        subprocess.run(
            [pdftoppm, "-png", "-r", str(PNG_DPI), "-singlefile",
             str(pdf_path), str(pdf_path.with_suffix(""))],
            check=True,
        )
        return png_path

    try:
        import pymupdf  # PyMuPDF >= 1.24
    except ModuleNotFoundError:
        try:
            import fitz as pymupdf  # older PyMuPDF
        except ModuleNotFoundError:
            return None

    with pymupdf.open(pdf_path) as doc:
        page = doc.load_page(0)
        page.get_pixmap(dpi=PNG_DPI).save(png_path)
    return png_path


def render() -> list[Path]:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    for basename, builder, title in DIAGRAMS:
        xdsm = builder()
        xdsm.write(basename, build=False, cleanup=True, outdir=str(OUTPUT_DIR))
        tex_path = OUTPUT_DIR / f"{basename}.tex"
        written.extend([OUTPUT_DIR / f"{basename}.tikz", tex_path])
        print(f"  {title}: wrote {basename}.tikz and {basename}.tex")

        pdf_path = _render_pdf(tex_path)
        if pdf_path is not None:
            written.append(pdf_path)
            png_path = _render_png(pdf_path)
            if png_path is not None:
                written.append(png_path)

    return written


def dump_matrix() -> str:
    """Text rendering of the DSM grid, for reviewing a draft without LaTeX.

    Laid out exactly as pyXDSM draws it: systems on the diagonal, and a
    coupling from src to target at (row = src, column = target).  That is the
    XDSM convention - data leaves a component HORIZONTALLY along its own row
    and enters the next one VERTICALLY down its column.

    So a mark ABOVE the diagonal is feed-forward and a mark BELOW it is
    feedback.  Counting the marks below the diagonal is the quickest check
    that the component ordering is right: there should be exactly as many as
    there are real loops, and no more.
    """
    lines: list[str] = []
    for basename, builder, title in DIAGRAMS:
        xdsm = builder()
        names = [s[0] for s in xdsm.systems]
        idx = {n: i for i, n in enumerate(names)}
        n = len(names)
        grid = [["  ." for _ in range(n)] for _ in range(n)]
        for i, name in enumerate(names):
            grid[i][i] = f"{i + 1:>3}"
        feedback = []
        for conn in xdsm.connections:
            src, tgt = conn[0], conn[1]
            r, c = idx[src], idx[tgt]
            if r < c:
                grid[r][c] = "  >"       # above the diagonal: feed-forward
            else:
                grid[r][c] = "  <"       # below the diagonal: feedback
                feedback.append((src, tgt))

        width = max(len(x) for x in names) + 1
        lines.append("=" * 78)
        lines.append(f"{title}   ({basename})")
        lines.append("=" * 78)
        header = " " * (width + 2) + "".join(f"{j + 1:>4}" for j in range(n))
        lines.append(header)
        for i, name in enumerate(names):
            lines.append(f"{name:>{width}}  " + "".join(grid[i]))
        lines.append("")
        lines.append(f"  systems {n}   couplings {len(xdsm.connections)}   "
                     f"feedback {len(feedback)}")
        for src, tgt in feedback:
            lines.append(f"    loop closure: {src} -> {tgt}")
        lines.append("")
    return "\n".join(lines)


def main() -> int:
    if "--matrix" in sys.argv:
        print(dump_matrix())
        return 0

    print("Building AAMF XDSM diagrams...")
    try:
        outputs = render()
    except MissingPrerequisiteError as exc:
        print(exc, file=sys.stderr)
        return 1

    matrix_path = OUTPUT_DIR / "xdsm_matrix_summary.txt"
    matrix_path.write_text(dump_matrix(), encoding="utf-8")
    outputs.append(matrix_path)
    print(f"  Matrix summary: wrote {matrix_path.name}")

    if not any(p.suffix == ".pdf" for p in outputs):
        print(
            "\nNo PDF produced: 'pdflatex' is not on PATH.\n"
            "The .tex and .tikz sources are complete and render anywhere:\n"
            "  winget install MiKTeX.MiKTeX      then reopen the shell and rerun\n"
            "or upload each .tex with its .tikz to Overleaf and compile there.\n"
            "PNG needs no extra tools once pdflatex exists: PyMuPDF is used if\n"
            "poppler's pdftoppm is absent (pip install PyMuPDF).",
            file=sys.stderr,
        )
    elif not any(p.suffix == ".png" for p in outputs):
        print(
            "\nPDF produced but no PNG: install either poppler (pdftoppm) or\n"
            "PyMuPDF (pip install PyMuPDF), then rerun.",
            file=sys.stderr,
        )
    print(f"\n{len(outputs)} files written to {OUTPUT_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

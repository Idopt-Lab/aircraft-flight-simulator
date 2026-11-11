# Silver Scythe — AAMF example case

Silver Scythe is the carrier strike-fighter concept of the Virginia Tech
Gobbler Ghost Program (AIAA final report, 2026): a 61,018 lb blended
lifting-body delta with 38° leading-edge sweep, twin F414-GE-400 engines
with two-axis thrust vectoring, and ruddervators in place of a
conventional tail.

This example puts that aircraft into the Aircraft-Agnostic Modular
Framework with two fidelity levels, built the same way as the F-16,
Boeing 737 and Cessna 172 cases: a data file, a coefficient lookup, a
builder, and a `*Simp.m` driver that trims the aircraft and publishes the
base-workspace variables `classes/FlightSim_Aircraft.slx` expects.

---

## Quick start

```matlab
cd AircraftDesign_Flightsim/MATLAB
addpath(genpath(pwd))

SilverScytheSelfCheck                     % verify every ported tool
out = SilverScytheSimp("high");           % trim + analyse + export
open_system(fullfile('classes','FlightSim_Aircraft.slx'))
```

Then run the Simulink model. For the CesiumJS view, start the web server
in `AircraftDesign_Flightsim/CesiumJS` (`python main.py`), open
`http://localhost:8000`, and run the model — the browser updates live.

`CesiumJS/main.js` renders whichever `.glb` it is pointed at; it ships
pointing at `f16.glb`, so until a Silver Scythe model is added the aircraft
flies as an F-16 shape. The trajectory, attitude and the whole
UDP → WebSocket data path are correct either way.

Other entry points:

```matlab
out  = SilverScytheSimp("low");           % textbook fidelity level
perf = SilverScythePerformance(SilverScytheContext("high"));
val  = SilverScytheValidation;            % low vs high vs report
```

---

## The two fidelity levels

The fidelity switch selects the **aerodynamic and propulsion methods**.
Everything else — geometry, mass build-up, frames, control surfaces,
thrust vectoring, the longitudinal derivative layer — is shared, so the
comparison isolates the modelling choice.

| | `low` | `high` |
|---|---|---|
| Section lift | thin-airfoil `2π`/rad, report `α_L0` = −1.05° | NACA 64-206 XFOIL polar, reduced by the same full-range fit as `cl_alpha.m` |
| 3-D lift slope | Nikolai/Raymer compressible swept finite wing (the `stall_solver.m` relation), Prandtl–Glauert with a slender-wing ceiling | same relation, XFOIL 2-D slope |
| Profile drag | `CD_min = Cfe·S_wet/S_ref`, Raymer Table 12.3 | **FRICTION** (Mason): per-component Eckert/Chapman–Rubesin laminar plus van Driest II turbulent friction with a Kármán–Schoenherr law, Raymer form factors, report Table 5 interference factors, 5% stealth leakage |
| Induced drag | Raymer swept Oswald efficiency, supersonic `k` correlation, 0.95–1.05 blend | **LIDRAG** (Ives/Mason): Trefftz-plane harmonic analysis of the 50-station spanload |
| Wave drag | Raymer `4.5π/S·(A_max/l)²·EWD` correlation | Raymer drag-divergence curve fits below Mach 1.2 plus the Sears-Haack estimate of report Eqs 15–16 |
| Engine | Mattingly algebraic deck, report §13.5 | **NPSS** deck, report §13.6, Mattingly-scaled outside its coverage |

The low level reproduces the senior-design textbook script
`drag_calc_sol.m`; the high level reproduces the actual report chain.

---

## Verification

`SilverScytheSelfCheck` compares each ported tool against a number the
original tool produced, not against another part of this package. All 30
checks pass:

| ported tool | model | source |
|---|---|---|
| FRICTION, `CD_profile` at Mach 0.80 / 30 kft | 0.0117085 | 0.0117085 (report Fig 24) |
| FRICTION, `CD_profile` at Mach 1.64 | 0.0092247 | 0.0092506 (−0.28%) |
| LIDRAG span efficiency | 0.88925 | 0.89127 (archived `lidrag.exe` run, −0.23%) |
| LIDRAG, elliptic spanload | 1.00000 | 1.0 exactly |
| Raymer wave drag, Mach 0.80 | 0.00051113 | 0.00051113 (report Fig 24) |
| Raymer wave drag, Mach 1.64 | 0.0161058 | 0.0161335 (−0.17%) |
| Total `CD` at Mach 1.64 | 0.0297839 | 0.0297827 (+0.004%) |
| NPSS design point | 13,700 lbf / 0.8552, 22,000 lbf / 1.8439 | report Table 11, exact |
| Neutral point | 22.719 ft | 22.7 ft (report Table 8) |
| Empty weight | 33,970 lb | report Table 12, exact |

The archived LIDRAG output is kept in `reference/lidrag_silverscythe.out`
so the span-efficiency check has an independent anchor.

---

## Performance against the report

`SilverScytheValidation` runs the low and high levels plus a
*report-reconstruction* case — the high-fidelity drag deck with the
Mattingly engine deck and the report's own lift bookkeeping, which is the
exact combination behind the published numbers.

| quantity | low | high | reconstruction | published |
|---|---:|---:|---:|---:|
| Max Mach, 30,000 ft, AB | 1.83 | 1.79 | **1.63** | 1.64 |
| Max Mach, sea level, MIL | 1.05 | 0.99 | 0.96 | 0.93 |
| Approach speed, kt | 137.3 | **134.4** | 134.4 | 133.8 |
| Arresting engagement speed, kt | 144.2 | **141.1** | 141.1 | 140 |
| Catapult end speed, kt | 136.8 | 136.8 | 136.8 | 133.5 |
| SEROC at launch, ft/min | 2304 | 2031 | **2186** | 2308 |
| SEROC at approach, ft/min | 5362 | 5001 | **3069** | 3320 |
| Sustained turn at 20,000 ft, °/s | 12.34 | **12.03** | 11.79 | 12 |
| Instantaneous turn at 20,000 ft, °/s | 21.60 | **21.66** | 21.65 | 22 |
| Corner speed at 20,000 ft, kt | 401.4 | **400.4** | 400.5 | 400 |
| Optimum nozzle deflection, ° | 11 | **11** | 11 | 11 |
| `CD` total, Mach 1.64 / 30 kft | 0.0220 | 0.0260 | **0.029784** | 0.029783 |

Turn performance, corner speed, the nozzle trade-study optimum and the
approach speeds land within 2% of the report on the high-fidelity model
directly. Maximum Mach and the rate of climb need the reconstruction case,
for the reason in finding 2 below.

---

## What the comparison turned up in the report

These are differences between the report and this implementation that are
**not** modelling choices — they are inconsistencies inside the source
data. Each is reported rather than tuned away.

**1 — Report Figure 24 is drawn at a CL about 27% above the level-flight
value.** `lift_calc.m` adds the induced-angle term of report Eq 13 to the
geometric angle instead of subtracting it, and multiplies by the 2-D
rather than the 3-D lift slope; because the normalising chord `c_a` is the
true average chord, the spanload integral then returns the section `cl`
as the aircraft `CL`. At the cruise point that is 0.345 against 0.272.
The error is confined to induced drag (0.01857 published against 0.01194
at the correct CL), and it is why the corrected model predicts Mach 1.79
rather than 1.64 and a higher single-engine climb rate at approach.
`SilverScytheReportCL` reinstates the original bookkeeping so the
published values can be reproduced; the flight model uses the CL from the
aircraft state, as it must.

**2 — Report Tables 8 and 9 disagree about how unstable the aircraft is.**
Table 8 gives a static margin of −0.23 MAC, and the ported
aerodynamic-centre chain confirms it from geometry (neutral point
22.719 ft against the published 22.7 ft). But at that margin the cruise
trim demand saturates the ±25° ruddervator limit of Table 7 — the
tail-volume pitch power of that empennage is only −0.10/rad, and thrust
vectoring adds about 0.005 in `Cm` at cruise thrust, so it cannot close
the gap. Table 9 (DATCOM), the higher-fidelity tool the team used to
validate the empirical chain, gives a near-neutral airframe that trims
with about 8° of ruddervator. The simulation uses the DATCOM layer;
`options.longitudinal_source = "textbook"` reproduces the Table 8 case and
its shortfall.

**3 — Report Table 13 does not sum to its own total.** The itemised
useful load adds to 27,018 lb against the 27,048 lb printed, and the
itemised stations give a takeoff CG of 31.074 ft against the printed
30.23 ft. The model follows the itemised entries, since those are the
traceable data, and reports the difference.

**4 — `emp_wave.m` does not match report Eq 15 as typeset.** The code
evaluates the Sears-Haack term as `(9π/2)(A_max/l)`, with the area ratio
to the first power, while Eq 15 squares it. The un-squared form is what
produced the published Figure 24 value, so it is the form implemented
here, with the discrepancy noted in `SilverScytheWaveDrag.m`.

**5 — The lateral-directional stability derivatives have no surviving
Silver Scythe high-fidelity source.** The DATCOM output files
(`SilverScytheAIL.out`, `SilverScytheFLAP.out`) were not archived with
the senior-design repository, and Table 9 read per radian leaves the
aircraft with `Cl_p` = −0.003 and `Cn_r` = +0.0002 — no roll or yaw
damping, and a destabilising yaw term. Control derivatives are computed
from geometry (exact in form, since the surfaces are all-moving), roll
damping from strip theory, and the remaining stability derivatives are
comparator-scaled from the F-16A deck in `examples/F16` by the ratio of
vertical tail-volume coefficients. Every term is tagged with its basis in
`ctx.sd.basis` and printed by `SilverScytheValidation`. To close them,
drop the original DATCOM export into `data/datcom_derivatives.csv` with
`Name,Value` columns.

---

## Control vector and thrust vectoring

```
1 aileron        wing trailing-edge differential   ±25°
2 elevator       symmetric ruddervator             ±25°
3 rudder         differential ruddervator          ±25°
4 nozzle_pitch   thrust-vector pitch command       ±20°
5 nozzle_yaw     thrust-vector yaw command         ±20°
6 throttle_L     left F414                          0–1
7 throttle_R     right F414                         0–1
```

Channels 2 and 3 keep the framework's generic names so the standard trim
configuration used by the other examples works unchanged; each surface
records its physical identity in `surface_type`.

Throttle convention matches the report's engine decks: 0 = idle,
0.5 = maximum dry (military), 1.0 = maximum afterburning.

Thrust vectoring carries **no aerodynamic derivative**. Channels 4 and 5
rotate each engine's thrust vector inside
`SilverScytheVectoringTurbofan`, and the force is applied at the nozzle
station (47.9 ft, report Table 7). That is the honest representation of
an effector whose moment scales with thrust rather than with dynamic
pressure, and it is why it keeps working at the low dynamic pressures of
carrier launch and recovery. Separate left and right throttles give the
differential thrust of report §13.1.

---

## Files

| file | role |
|---|---|
| `SilverScytheData.m` | every input, each tagged `[RPT]` report, `[SD]` senior-design deck, or `[SIM]` simulator-only |
| `SilverScytheAirfoil.m` | NACA 64-206 XFOIL polar reader, port of `cl_alpha.m` |
| `SilverScytheFriction.m` | FRICTION (Mason) profile and form drag |
| `SilverScytheLidrag.m` | LIDRAG (Ives/Mason) Trefftz-plane induced drag |
| `SilverScytheSpanload.m` | 50-station spanload, port of `lift_calc.m` |
| `SilverScytheWaveDrag.m` | Raymer/Brandt wave drag, port of `emp_wave.m` |
| `SilverScytheAeroDeck.m` | builds and caches the high-fidelity deck |
| `SilverScytheLift.m` | lift curve for both fidelities |
| `SilverScytheDrag.m` | drag polar for both fidelities |
| `SilverScytheStabDerivatives.m` | S&C derivative layer, every term tagged with its basis |
| `SilverScytheEngineDeck.m` | NPSS deck loader |
| `SilverScytheThrustModel.m` | Mattingly and NPSS thrust/TSFC |
| `SilverScytheVectoringTurbofan.m` | thrust-vectoring propulsive element |
| `SilverScytheContext.m` | assembles one fidelity level |
| `SilverScytheLookup.m` | AAMF coefficient callback |
| `SilverScytheBuild.m` | builds the `Aircraft` object |
| `SilverScytheSimp.m` | **main driver**: trim, analyses, Simulink export |
| `SilverScythePerformance.m` | report §16 performance analyses |
| `SilverScytheReportCL.m` | the report's own CL bookkeeping, for reproduction only |
| `SilverScytheSelfCheck.m` | 30 checks of the ported tools against their sources |
| `SilverScytheValidation.m` | low vs high vs report, writes `results/` |
| `data/NACA_64206_xfoil.txt` | XFOIL polar, Re 12e6 |
| `data/F414_NPSS_deck.csv` | 1,039 NPSS cases from `F414-GE-400_Outfile.txt` |
| `reference/lidrag_silverscythe.out` | archived `lidrag.exe` output |
| `results/` | generated aero deck cache and validation output |

The aero deck is cached in `results/SilverScytheAeroDeck.mat` and rebuilt
automatically whenever the geometry or grid signature changes; force a
rebuild with `SilverScytheContext("high",struct('rebuild_deck',true))`.

---

## Sources

* AIAA final report, *Gobbler Ghost Program*, Virginia Tech, 2026 —
  geometry, weights, control surfaces, published performance
* `Aero/SUMAC/` — `zero_lift_drag.m` (FRICTION), `induced_drag.m` and
  `lift_calc.m` (LIDRAG), `emp_wave.m`, `SUMAC.m`
* `Performance/` — `drag_calc_sol.m`, `stall_solver.m`,
  `takeoff_landing.m`, `turn_performance.m`, `climb_func.m`,
  `FlightPerformance.m`, `breguet.m`
* `Propulsion/` — `computeThrust414.m`, `F414-GE-400_Outfile.txt` (NPSS)
* `S&C/SMcalcs/` — `AeroCenterChord*.m`, `AerodynamicCenter.m`,
  `NeutralPointStaticMargin.m`
* `Req_Analysis/aircraft_params.m`, `W&B/calculateWeightBalance.m`
* Mason, W. H., FRICTION; Ives, D. and Mason, W. H., LIDRAG
* Raymer, *Aircraft Design: A Conceptual Approach*; Nicolai and Carichner;
  Roskam, *Airplane Design*; Mattingly, *Aircraft Engine Design*

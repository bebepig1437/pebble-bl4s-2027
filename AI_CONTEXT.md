# AI_CONTEXT.md — CERN BL4S 2027 Bismuth Scintillator Simulation

> **Instructions for AI Assistants (IDE / Browser):** Read this file completely before modifying any code or writing proposal text. It contains verified 2027 competition constraints, corrected physics assumptions, and safe-editing rules.

## 1. Project Identity & Scientific Goal
- **Competition:** CERN Beamline for Schools (BL4S) 2027 (Proposals due March 2027).
- **Theme Alignment:** "Build and test your own detector" & "Study different ways to detect photons."
- **Goal:** Benchmark lead-free, low-cost (<$3.20/disc) bismuth-based composite scintillators against commercial plastic (`EJ200`) under relativistic single-electron beams, characterizing light yield (photoelectrons/MIP), timing kinetics, multiple scattering, and Bremsstrahlung output.
- **Candidate Formulations:**
  1. `BiNanoPDMS_25wt` (Method B, realistic wet-cast): 25 wt% Sb-doped 3D double perovskite nanocrystals ($Cs_2AgBi_{0.9}Sb_{0.1}Br_6$, $\rho_c = 4.80\text{ g/cm}^3$, $d < 20\text{ nm}$ Rayleigh regime) in PDMS Sylgard 184 ($C_2H_6OSi$, $\rho_m = 1.02\text{ g/cm}^3$). Exact rule-of-mixtures density = **1.27 g/cm³**.
  2. `BiNanoPDMS_66wt` (Method B, high-load target): 66 wt% (~29 vol%) loading to match the ~2.12 g/cm³ composite density cited in early project dossiers.
  3. `BiPh3PVK` (Method C, fast benchmark): 15 wt% Triphenylbismuth ($C_{18}H_{15}Bi$) dissolved molecularly in PVK ($C_{14}H_{11}N$) with FRET fluors ($\tau_1 \approx 1.6\text{ ns}$, $\rho \approx 1.25\text{–}1.30\text{ g/cm}^3$).
  4. `EJ200` (Control): Commercial polyvinyltoluene plastic scintillator ($\rho = 1.023\text{ g/cm}^3$, $\tau = 2.1\text{ ns}$, $Y = 10,000\text{ ph/MeV}$).
  5. `LYSO` (Optional reference): $Lu_{1.8}Y_{0.2}SiO_5$ ($\rho = 7.10\text{ g/cm}^3$, intrinsic $^{176}Lu$ noise ~40 Bq/g, >$300/cm³).

## 2. MANDATORY 2027 Facility Constraints (DO NOT VIOLATE)
- **NO CERN PS / NO Hadrons / NO Muons / NO DWCs:** In 2027–2028, the CERN PS accelerator is shut down for maintenance. Only electron and positron beams are available across CERN, DESY, and ELSA.
- **DO NOT Target CERN CLEAR for Single-Particle Timing:** CLEAR delivers 60–220 MeV/c dense electron bunches ($\ge 6.24 \times 10^7$ electrons per 0.01 nC bunch), making single-particle Landau/CTR measurements impossible.
- **Primary Target Facility — DESY-II (Hamburg):**
  - Pure `e-` or `e+` single-particle beam, `1–6 GeV/c`, momentum spread `~0.15 GeV/c`, round beam spot `20 mm` diameter.
  - Available detectors: Scintillator paddles, 6-plane Silicon Beam Telescope ($2\times 2\text{ cm}^2$, few $\mu\text{m}$ resolution), 16 Lead-Glass Calorimeters ($10\times 10\times 37\text{ cm}$, $\sigma_E/E = 0.02\% + 6.3\%/\sqrt{E}$), Big Red Dipole Magnet (1.35 T).
- **Secondary Target Facility — ELSA (University of Bonn):**
  - Pure `e-` beam, `0.8–3.2 GeV` (2.9 or 3.2 GeV in secondary user mode), ANEMONE 6-plane Mimosa26 silicon telescope ($18.4\ \mu\text{m}$ pitch), ATLAS FE-I4 timing plane, scintillators.

## 3. Corrected Physics & Safety Rules (Guardrails Against Old Dossier Errors)
1. **No $Z^4$ Photoelectric Claims for GeV Electrons:** Photoelectric absorption ($\propto \rho Z^4/E^3$) governs keV X-rays, NOT 1–6 GeV electrons. At GeV energies, energy loss is governed by **Bethe-Bloch ionization** (scales with $\rho \cdot Z/A$), **Multiple Coulomb Scattering** (Highland formula $\propto \sqrt{x/X_0}$), and **Bremsstrahlung** ($X_0^{-1} \propto \rho Z^2/A$).
2. **No "25 ns LHC Bunch Pileup" Claim at DESY-II:** DESY-II delivers single electrons continuously over an 80 ms cycle. Frame timing around single-electron **Coincidence Time Resolution (CTR)** and fast scintillation decay kinetics ($\tau \approx 15\text{–}30\text{ ns}$ via 10 mol% $Sb^{3+}$ symmetry breaking vs. >1,000 ns undoped STEs).
3. **Rule-of-Mixtures Density Honesty:** Never hardcode $\rho = 2.12\text{ g/cm}^3$ for 25 wt% loading. Always compute $\rho = (w_c/\rho_c + (1-w_c)/\rho_m)^{-1}$ dynamically.
4. **Chemistry Safety:** Hot 48 wt% HBr (110–120 °C) and $Sb_2O_3$ MUST be synthesized in a certified fume hood with a university partner lab—never claim "no fume hood needed."

## 4. Hardware Parameters (Readout & Optics)
- **DUT Disc:** 25 mm diameter, 1–3 mm thick (nominal 2 mm).
- **Trigger Paddles (`S1`, `S2`):** EJ-200, $10\times 10\times 2\text{ mm}$ at $z = -60\text{ mm}$ and $z = +60\text{ mm}$. Hardware coincidence cut: `edep_S1 > 0.05 MeV && edep_S2 > 0.05 MeV`.
- **SiPM:** Hamamatsu S13360-3050CS ($3.0\times 3.0\text{ mm}^2$, $50\ \mu\text{m}$ pitch, 3,600 pixels, $V_{br} \approx 53\text{ V}$, +3.0 V overvoltage, **PDE = 40% at 450 nm**).
- **Optical Interface:** Dow Corning Q2-3067 grease ($n = 1.465$), SiPM epoxy ($n = 1.55$), 3M ESR specular reflector foil ($R = 0.985$ on all non-SiPM faces).
- **Light Yield Estimates:** EJ-200 = `10,000 ph/MeV`; Sb-doped $Cs_2AgBiBr_6$ = `5,000–15,000 ph/MeV` (nominal `10,000`); $BiPh_3$-PVK = `6,000–8,000 ph/MeV` (nominal `7,000`).

## 5. Codebase Architecture & Safe-Editing Rules
- **`config.json`**: Edit beam energies, materials, thicknesses, event counts, and optical properties HERE. Do not hardcode numbers in `run_pipeline.py` or `src/sim.cc`.
- **`src/sim.cc`**: Single-file serial Geant4 executable using `FTFP_BERT` + `G4EmStandardPhysics_option4`. Controlled strictly via CLI flags (`--material`, `--loading`, `--thickness`, `--energy`, `--espread`, `--particle`, `--events`, `--seed`, `--out`, `--vis`). Outputs clean CSV + `.meta.json` (containing Geant4's exact calculated density and radiation length $X_0$).
- **`src/optical_mc.py`**: Vectorized 3D NumPy photon ray-tracer. Handles Fresnel/TIR reflection, 3M ESR walls, bulk absorption, and Rayleigh vs. Mie scattering. **Do NOT move optical ray-tracing back into Geant4 C++** (`G4OpticalPhysics` is 100x slower and fragile across Geant4 versions).
- **`run_pipeline.py`**: Compiles `src/sim.cc`, runs Geant4 sweeps in parallel via `ProcessPoolExecutor`, calls `optical_mc.py`, fits Moyal/Landau distributions, validates against Highland MCS theory, and writes `results/fig1..fig4` + `results/summary_table.md`.
- **Testing Rule:** Always verify changes with `python run_pipeline.py --smoke-test` (fast 2,000-event run) before running full statistics.

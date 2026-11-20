# CERN BL4S 2027 — Lead-Free Bismuth Scintillator Simulation Suite

**PEBBLE** — **P**erovskite **E**lectron **B**eam **B**enchmark for **L**ight **E**mission — is the name of the project and of this repository (`pebble-bl4s-2027`). Six letters, six words: *for* is the only word that contributes none. Each letter names a load-bearing piece of the simulation:

| Letter | Word | What it names here |
| --- | --- | --- |
| **P** | Perovskite | the lead-free $Cs_2AgBi_{0.9}Sb_{0.1}Br_6$ double perovskite dispersed in PDMS, the material the beam test is built around |
| **E** | Electron | the single-electron DESY-II (1–6 GeV/c) and ELSA (0.8–3.2 GeV) test beams |
| **B** | Beam | the 20 mm beam spot and $S_1/S_2$ paddle geometry this code reproduces (≈31% trigger acceptance) |
| **B** | Benchmark | the reference numbers — MPV $\Delta E$, $\epsilon_{\text{geo}}$, $N_{\text{p.e.}}$, $\theta_{\text{RMS}}$ — the experiment has to reproduce |
| **L** | Light | scintillation photons traced through total internal reflection, Rayleigh/Mie scattering and 3M ESR foil reflection |
| **E** | Emission | the photoelectron yield at the SiPM, the number that decides whether a clean trigger is possible |

Reproducible **Geant4 + 3D Optical Monte Carlo** digital twin for our **CERN Beamline for Schools (BL4S) 2027** proposal, configured for the **DESY-II (1–6 GeV/c)** and **ELSA (0.8–3.2 GeV)** single-electron test beamlines.

## What This Simulation Proves
1. **Energy Deposition (Landau-Vavilov MPV):** Exact thin-target relativistic electron energy loss ($\Delta E$) across `EJ-200`, `BiNanoPDMS` ($Cs_2AgBi_{0.9}Sb_{0.1}Br_6$ in PDMS), and `BiPh3-PVK`.
2. **Optical Signal Viability (SiPM Photoelectrons):** 3D Monte Carlo photon ray-tracing through a 25 mm disc wrapped in 3M ESR reflector foil (98.5%) into a $3\times 3\text{ mm}^2$ Hamamatsu S13360-3050CS SiPM (40% PDE), proving we clear the **>20 p.e.** clean trigger threshold and quantifying the **Rayleigh (Method B nanocrystals) vs. Mie (Method A microcrystals)** light-collection advantage.
3. **High-$Z$ Beamline Signatures:** Multiple Coulomb Scattering angles ($\theta_{\text{RMS}}$, cross-checked against the analytical Highland formula) and Bremsstrahlung photon production measurable by DESY-II's 6-plane silicon beam telescope and lead-glass calorimeter.

---

## Quickstart (Mac, Linux, or Windows WSL2)

### 1. Create the Reproducible Environment (One-Time Setup)
```bash
conda env create -f environment.yml
conda activate bl4s-sim
```

> **Environment note.** `environment.yml` pins the headless `geant4=*=noqt*` build — no Qt or
> OpenGL viewer, so `--vis` reports that it is unavailable (see step 5). It also lists the
> `expat` and `zlib` dev packages: Geant4's CMake config runs `find_dependency(EXPAT 2.8.5)`
> and `(ZLIB 1.3.2)`, and on macOS CMake otherwise selects the system SDK stubs, whose older
> versions fail that version check. The interactive viewer needs a `qt` Geant4 build plus
> `cmake -DBL4S_VIS=ON`.

### 2. Run a Fast Verification Smoke Test (~30 seconds)
```bash
python run_pipeline.py --smoke-test
```

### 3. Run the Full High-Statistics Simulation & Generate All Plots (~2–3 minutes)
```bash
python run_pipeline.py
```

### 4. Re-Plot or Tweak Optical/SiPM Settings Without Re-Running Geant4 (~3 seconds)
```bash
python run_pipeline.py --skip-geant4
```

### 5. Optional: Launch Interactive 3D Geant4 Viewer
```bash
./build/bl4s_sim --material BiNanoPDMS --loading 0.25 --thickness 2.0 --energy 3.0 --vis
```

---

## Outputs Generated in `results/`

| File | Description |
| --- | --- |
| `fig1_setup_3d.png` | 3D engineering render of the DESY-II telescope, $S_1/S_2$ triggers, DUT disc, SiPM, and calorimeter |
| `fig2_landau_edep.png` | Overlaid 3 GeV electron energy-deposition spectra (2 mm disc) with Moyal/Landau MPV fits |
| `fig3_sipm_yield.png` | Expected Hamamatsu SiPM photoelectrons vs. thickness (with literature light-yield uncertainty bands & 20 p.e. trigger floor) + Rayleigh vs. Mie collection inset |
| `fig4_scattering_brem.png` | RMS electron scattering angle vs. analytical Highland theory & Bremsstrahlung photon yield across 1–6 GeV |
| `summary_table.md` | Proposal-ready Markdown benchmark table with MPV $\Delta E$, $X_0$, $\epsilon_{\text{geo}}$, $N_{\text{p.e.}}$, and $\theta_{\text{RMS}}$ |

---

## How to Edit Parameters Safely

Open **`config.json`** to change:

* Beam facility (`DESY-II` or `ELSA`), particle (`e-` or `e+`), energies (`1, 3, 6 GeV`), or momentum spread
* DUT disc diameter, thickness sweep (`1, 2, 3 mm`), and composite weight-fraction loading
* Optical constants (refractive index, absorption/scattering lengths, 3M ESR reflectivity, SiPM PDE)

No C++ editing is required for standard parameter sweeps.

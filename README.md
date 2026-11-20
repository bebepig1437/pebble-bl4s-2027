# PEBBLE — Lead-Free Bismuth Scintillators for CERN BL4S 2027

PEBBLE (**P**erovskite **E**lectron **B**eam **B**enchmark for **L**ight **E**mission) is a student project for the [CERN Beamline for Schools](https://beamline-for-schools.web.cern.ch/) competition. It asks one question: on a GeV electron beam, can a lead-free bismuth-loaded plastic scintillator beat commercial EJ-200 in a detector only a few millimetres thick?

This repository is the simulation work behind that question — a Geant4 model of energy deposition and a ray-traced model of light collection, both driven from one config file.

## The question

A 2 mm EJ-200 disc sits at 0.47% of a radiation length. It is fast and cheap, but thin and light, so a minimum-ionising 3 GeV electron deposits only about 0.3 MeV in it.

Loading the polymer with bismuth changes that. The composite gets denser and its radiation length shortens, so each millimetre absorbs more energy and produces more bremsstrahlung and scattering. But the optics get worse as the loading goes up: the material becomes hazier, less light reaches the SiPM, and the molecular BiPh3-PVK route starts from a lower light yield.

So the interesting question is not whether bismuth is "better", but where the product of deposited energy and collected light peaks — and whether any of these materials clears the photoelectron threshold needed for a clean trigger.

## The setup

![Simulated beamline](results/fig1_setup_3d.png)

| Component | Description |
| --- | --- |
| Beam | DESY-II single electrons, 1–6 GeV/c, 20 mm spot; ELSA (0.8–3.2 GeV) as the alternative |
| Trigger | Two 10 × 10 × 2 mm EJ-200 paddles at z = ∓60 mm, coincidence on > 0.05 MeV in each |
| Sample | 25 mm disc, 1–3 mm thick (2 mm nominal), wrapped in 3M ESR foil (R = 98.5%) |
| Readout | Hamamatsu S13360-3050CS SiPM, 3 × 3 mm², 40% PDE at 450 nm, coupled with n = 1.465 grease |
| Threshold | 20 photoelectrons for a clean trigger |

The 20 mm beam and the 10 mm paddles give roughly 31% trigger acceptance, which matches the real DESY-II beam profile rather than an idealised pencil beam.

## What the simulation does

**Energy deposition (Geant4).** A single-file C++ application built on `FTFP_BERT` with `G4EmStandardPhysics_option4` tracks each electron through the paddles, the sample and the calorimeter, recording energy deposited in every volume. Composite densities are computed from the rule of mixtures and cross-checked against Geant4's own material calculation, never hardcoded.

**Light collection (NumPy ray tracer).** Scintillation photons are traced from their production point through Fresnel and total-internal reflection at the disc faces, absorption and Rayleigh or Mie scattering in the bulk, and specular reflection off the ESR wrap, until they hit the SiPM or are lost. This runs in vectorised NumPy rather than `G4OpticalPhysics`, which keeps a 100,000-photon sample to a fraction of a second.

**Checks against theory.** Energy-loss spectra are fitted with a Landau/Moyal most-probable value, and the simulated multiple-scattering angle is compared with the analytical Highland formula. Where the two disagree, the disagreement is reported rather than smoothed over.

## Results

3 GeV electrons on a 2 mm disc, 25,000 events per material, DESY-II geometry:

| Material | Density (g/cm³) | X₀ (mm) | ΔE MPV (MeV) | Light collection | SiPM yield (p.e.) |
| --- | ---: | ---: | ---: | ---: | ---: |
| EJ-200 (control) | 1.03 | 425 | 0.315 | 19.0% | 239 |
| Bi-nano PDMS, 25 wt% | 1.27 | 151 | 0.362 | 22.9% | 332 |
| Bi-nano PDMS, 66 wt% | 2.12 | 56 | 0.552 | 17.7% | 391 |
| BiPh3-PVK, 15 wt% | 1.25 | 245 | 0.362 | 15.2% | 154 |

![Photoelectron yield versus thickness](results/fig3_sipm_yield.png)

Three things follow from this table.

The 66 wt% composite deposits 1.75× the energy of EJ-200 and gives the highest yield at 391 p.e., but its collection efficiency is the second worst of the four. The gain comes from density and radiation length, not from better optics — heavier loading makes the material hazier, exactly as expected.

At 25 wt%, the realistic wet-cast loading, the composite still returns 39% more photoelectrons than EJ-200 while keeping the best collection efficiency in the set. That is the formulation the proposal is built around.

BiPh3-PVK is the dimmest at 154 p.e. despite matching the 25 wt% composite in deposited energy, because its light yield per MeV is lower and its higher refractive index traps more light. All four clear the 20 p.e. threshold, so the choice is about margin rather than feasibility.

The optical model also separates the two synthesis routes: transparent Rayleigh-sized nanocrystals collect 22.9% of the light against 16.9% for hazy Mie-sized microcrystals at the same loading, a 1.4× difference that is purely optical.

## Running it

```bash
conda env create -f environment.yml
conda activate bl4s-sim

python run_pipeline.py --smoke-test    # 2,000 events, about 30 s
python run_pipeline.py                 # full statistics, about 2–3 min
python run_pipeline.py --skip-geant4   # re-plot from existing CSVs, about 3 s
```

`environment.yml` pins the headless `geant4=*=noqt*` build, so `--vis` reports that no viewer is available; a Qt build plus `-DBL4S_VIS=ON` enables it. The `expat` and `zlib` entries are there because Geant4's CMake config requires versions newer than the macOS SDK ships.

## Repository layout

| Path | Contents |
| --- | --- |
| `config.json` | Beam, geometry, material and optical parameters — the only file that needs editing for a sweep |
| `src/sim.cc` | Geant4 application: geometry, stepping action, CSV output |
| `src/optical_mc.py` | Vectorised photon transport and collection efficiency |
| `run_pipeline.py` | Builds the simulation, runs the sweep in parallel, fits, cross-checks and plots |
| `results/` | Figures, per-job metrics and the summary table |

## Scope and limits

This is a student project, not a measurement. Every number here is a simulation estimate intended to size a beam-test proposal, and the figures should be read that way.

The material parameters are literature values or informed estimates: light yields are taken as ranges and plotted with uncertainty bands, and the absorption and scattering lengths are nominal. The optical model is geometric and monochromatic at the 450 nm design point. Detector electronics, timing resolution and calibration are not simulated. The Highland comparison agrees to a few percent at low loading and departs at 66 wt%, where the Gaussian approximation is no longer a good description of a thick, high-Z sample.

What the simulation does support is a comparison: under one consistent set of assumptions, the bismuth-loaded discs return more photoelectrons than EJ-200, and the 25 wt% formulation is the best compromise between energy deposition and light collection.

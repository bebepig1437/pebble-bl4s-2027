#!/usr/bin/env python3
"""
Master orchestrator for the CERN BL4S 2027 bismuth scintillator simulation.

1. Builds src/sim.cc with CMake (incremental).
2. Runs the Geant4 sweeps from config.json in parallel.
3. Runs the vectorized 3D optical photon Monte Carlo from src/optical_mc.py.
4. Fits Moyal/Landau MPVs, cross-checks Geant4 scattering against the Highland
   formula, and writes results/fig1..fig4.png plus results/summary_table.md.
"""

import argparse
import json
import subprocess
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.stats import moyal

from src.optical_mc import simulate_optical_collection

ROOT = Path(__file__).resolve().parent
BUILD_DIR = ROOT / "build"
RESULTS_DIR = ROOT / "results"
DATA_DIR = RESULTS_DIR / "raw_csv"
BIN_PATH = BUILD_DIR / "bl4s_sim"


# ---------------------------------------------------------------------------
# Sweep definition
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Beam:
    """Beamline constants shared by every job."""

    particle: str
    spread_GeV: float
    spot_mm: float
    disc_mm: float


@dataclass(frozen=True)
class Job:
    """One Geant4 run: a material sample in the beam."""

    mat_id: str
    material: str
    loading: float
    thickness: float
    energy: float
    events: int
    seed: int
    out_csv: Path
    beam: Beam

    @property
    def tag(self) -> str:
        return f"{self.mat_id}_t{self.thickness:.1f}mm_E{self.energy:.1f}GeV"


def build_jobs(cfg, events):
    """Thickness scan at nominal energy, then energy scan at nominal thickness."""
    beam = Beam(
        particle=cfg["project"]["particle"],
        spread_GeV=cfg["beam"]["energy_spread_GeV"],
        spot_mm=cfg["beam"]["beam_diameter_mm"],
        disc_mm=cfg["geometry"]["disc_diameter_mm"],
    )
    nom_E = cfg["beam"]["nominal_energy_GeV"]
    nom_T = cfg["geometry"]["nominal_thickness_mm"]

    jobs, seen = [], set()
    seed = cfg["project"]["seed"]
    for mat in cfg["materials"]:
        samples = [(t, nom_E) for t in cfg["geometry"]["thickness_scan_mm"]]
        samples += [(nom_T, e) for e in cfg["beam"]["energy_scan_GeV"]]
        for thickness, energy in samples:
            key = (mat["id"], thickness, energy)
            if key in seen:
                continue
            seen.add(key)
            tag = f"{mat['id']}_t{thickness:.1f}mm_E{energy:.1f}GeV"
            jobs.append(Job(
                mat_id=mat["id"],
                material=mat["sim_material"],
                loading=mat["loading_wt"],
                thickness=thickness,
                energy=energy,
                events=events,
                seed=seed,
                out_csv=DATA_DIR / f"{tag}.csv",
                beam=beam,
            ))
            seed += 17
    return jobs


# ---------------------------------------------------------------------------
# Geant4
# ---------------------------------------------------------------------------
def build():
    """Configure and build the Geant4 executable (CMake is incremental)."""
    print("[1/4] Building the Geant4 executable...")
    subprocess.run(["cmake", "-S", str(ROOT), "-B", str(BUILD_DIR)], check=True)
    subprocess.run(["cmake", "--build", str(BUILD_DIR), "--parallel"], check=True)


def run(job):
    subprocess.run(
        [
            str(BIN_PATH),
            "--material", job.material,
            "--loading", f"{job.loading:.4f}",
            "--thickness", f"{job.thickness:.3f}",
            "--diameter", f"{job.beam.disc_mm:.3f}",
            "--energy", f"{job.energy:.3f}",
            "--espread", f"{job.beam.spread_GeV:.3f}",
            "--beam-diam", f"{job.beam.spot_mm:.3f}",
            "--particle", job.beam.particle,
            "--events", str(job.events),
            "--seed", str(job.seed),
            "--out", str(job.out_csv),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    return job.tag, job.out_csv


# ---------------------------------------------------------------------------
# Optics and analysis helpers
# ---------------------------------------------------------------------------
def optical_efficiency(cfg, thickness, rindex, abs_mm, scat_mm, n_photons, seed):
    o, g = cfg["optical_mc"], cfg["geometry"]
    return simulate_optical_collection(
        disc_diameter_mm=g["disc_diameter_mm"],
        thickness_mm=thickness,
        trigger_xy_mm=g["trigger_xy_mm"],
        sipm_xy_mm=g["sipm_xy_mm"],
        rindex=rindex,
        abs_length_mm=abs_mm,
        scat_length_mm=scat_mm,
        grease_rindex=o["grease_rindex"],
        sipm_epoxy_rindex=o["sipm_epoxy_rindex"],
        esr_reflectivity=o["esr_reflectivity"],
        n_photons=n_photons,
        max_bounces=o["max_bounces"],
        seed=seed,
    )


def highland_mrad(energy_GeV, x_over_X0):
    """PDG Highland formula for the planar RMS scattering angle [mrad]."""
    if x_over_X0 <= 0:
        return 0.0
    theta0_rad = (13.6 / (energy_GeV * 1000.0)) * np.sqrt(x_over_X0) * (1.0 + 0.038 * np.log(x_over_X0))
    return float(theta0_rad * 1000.0)


def fit_moyal_mpv(edep_MeV):
    """Fit the Landau-approximating Moyal distribution; its loc parameter is the mode."""
    valid = edep_MeV[(edep_MeV > 0.01) & (edep_MeV < np.percentile(edep_MeV, 98.5))]
    loc, scale = moyal.fit(valid)
    return float(loc), float(scale), float(scale / np.sqrt(max(len(valid), 1)))


def collect_metrics(cfg, jobs, optics, nom_E, nom_T):
    """Per-job physics metrics plus the nominal spectra used by the figures."""
    opt_cfg = cfg["optical_mc"]
    pde = opt_cfg["sipm_pde_450nm"]
    materials = {m["id"]: m for m in cfg["materials"]}

    rows, nominal = [], {}
    for job in jobs:
        if not job.out_csv.exists():
            raise SystemExit(f"Missing {job.out_csv}; run without --skip-geant4 first.")
        meta = json.loads(Path(f"{job.out_csv}.meta.json").read_text())
        trig = pd.read_csv(job.out_csv).query("triggered == 1")
        mpv, scale, mpv_err = fit_moyal_mpv(trig["edep_DUT_MeV"].to_numpy())

        # Planar RMS, trimmed at the 98% core to suppress hard Rutherford tails
        # (the regime the Highland formula describes).
        sx = trig["scatter_x_mrad"].to_numpy()
        sy = trig["scatter_y_mrad"].to_numpy()
        core_x = sx[np.abs(sx) < np.percentile(np.abs(sx), 98.0)]
        core_y = sy[np.abs(sy) < np.percentile(np.abs(sy), 98.0)]

        mat = materials[job.mat_id]
        eps_geo = optics[(job.mat_id, job.thickness)]
        rows.append({
            "mat_id": job.mat_id,
            "label": mat["label"],
            "color": mat["color"],
            "thickness_mm": job.thickness,
            "energy_GeV": job.energy,
            "n_triggered": len(trig),
            "density_g_cm3": meta["density_g_cm3"],
            "radlen_mm": meta["radlen_mm"],
            "x_over_X0_pct": 100.0 * meta["x_over_X0"],
            "mpv_MeV": mpv,
            "moyal_scale": scale,
            "mpv_err_MeV": mpv_err,
            "eps_geo_pct": 100.0 * eps_geo,
            "n_pe_nom": mpv * mat["light_yield_nominal"] * eps_geo * pde,
            "n_pe_min": mpv * mat["light_yield_min"] * eps_geo * pde,
            "n_pe_max": mpv * mat["light_yield_max"] * eps_geo * pde,
            "rms_scatter_mrad": float(np.sqrt(0.5 * (np.mean(core_x**2) + np.mean(core_y**2)))),
            "highland_mrad": highland_mrad(job.energy, meta["x_over_X0"]),
            "brem_rate_pct": float(100.0 * np.mean(trig["brem_photons"] > 0)),
        })
        if job.thickness == nom_T and job.energy == nom_E:
            nominal[job.mat_id] = trig

    return pd.DataFrame(rows), nominal


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def plot_setup_3d(cfg):
    """3D schematic of the DESY-II telescope, triggers, DUT, SiPM and calorimeter."""
    fig = plt.figure(figsize=(10, 5.2), dpi=300)
    ax = fig.add_subplot(111, projection="3d")

    def draw_box(z_center, dx, dy, dz, color, alpha=0.55):
        x0, x1 = -dx / 2, dx / 2
        y0, y1 = -dy / 2, dy / 2
        z0, z1 = z_center - dz / 2, z_center + dz / 2
        verts = [
            [(z0, x0, y0), (z0, x1, y0), (z0, x1, y1), (z0, x0, y1)],
            [(z1, x0, y0), (z1, x1, y0), (z1, x1, y1), (z1, x0, y1)],
            [(z0, x0, y0), (z1, x0, y0), (z1, x1, y0), (z0, x1, y0)],
            [(z0, x0, y1), (z1, x0, y1), (z1, x1, y1), (z0, x1, y1)],
            [(z0, x0, y0), (z1, x0, y0), (z1, x0, y1), (z0, x0, y1)],
            [(z0, x1, y0), (z1, x1, y0), (z1, x1, y1), (z0, x1, y1)],
        ]
        ax.add_collection3d(Poly3DCollection(verts, facecolors=color, edgecolors="k", linewidths=0.4, alpha=alpha))

    draw_box(-60, 10, 10, 2, "#2563EB", 0.7)   # S1
    draw_box(+60, 10, 10, 2, "#2563EB", 0.7)   # S2
    for zt in [-42, -30, -18, 18, 30, 42]:     # silicon beam telescope
        draw_box(zt, 20, 20, 0.8, "#64748B", 0.30)

    R = cfg["geometry"]["disc_diameter_mm"] / 2.0
    theta = np.linspace(0, 2 * np.pi, 48)
    xc, yc = R * np.cos(theta), R * np.sin(theta)
    for zface in (-1.0, +1.0):                 # DUT disc faces
        verts = [list(zip(np.full_like(xc, zface), xc, yc))]
        ax.add_collection3d(Poly3DCollection(verts, facecolors="#EA580C", edgecolors="#7C2D12", alpha=0.85))

    draw_box(1.5, 3, 3, 1.0, "#10B981", 0.95)  # SiPM
    draw_box(105, 35, 35, 30, "#CA8A04", 0.35)  # calorimeter front section

    rng = np.random.default_rng(7)
    for _ in range(12):                        # sample beam tracks
        x0, y0 = rng.uniform(-4, 4, 2)
        sx, sy = rng.normal(0, 0.02, 2)
        ax.plot([-90, 0, 105], [x0, x0, x0 + sx * 105], [y0, y0, y0 + sy * 105],
                color="#EF4444", linewidth=0.9, alpha=0.75)

    # Labels sit in two clearly separated bands (above / below the beam axis) so
    # the oblique projection cannot stack them on top of each other.
    above = dict(ha="center", va="bottom", fontsize=8)
    below = dict(ha="center", va="top", fontsize=8)
    ax.text(-74, 0, 22, "S1 Trigger\n(EJ-200)", color="#1D4ED8", weight="bold", **above)
    ax.text(58, 0, 22, "S2 Trigger", color="#1D4ED8", weight="bold", **above)
    ax.text(-80, 0, -22, "3x Si Telescope", color="#334155", **below)
    ax.text(8, 0, -22, "DUT Disc (25 mm)\n+ 3x3 mm SiPM", color="#C2410C", weight="bold", **below)
    ax.text(100, 0, -22, "Pb-Glass\nCalorimeter", color="#854D0E", weight="bold", **below)

    ax.set_xlim(-95, 125)
    ax.set_ylim(-25, 25)
    ax.set_zlim(-25, 25)
    ax.set_xlabel("Beam Axis Z [mm]")
    ax.set_ylabel("X [mm]")
    ax.set_zlabel("Y [mm]")
    ax.view_init(elev=22, azim=-62)
    ax.set_title(f"CERN BL4S 2027 ({cfg['project']['facility']}) Test-Beam Geometry — "
                 f"{cfg['beam']['nominal_energy_GeV']:.0f} GeV e⁻", pad=12, weight="bold")
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "fig1_setup_3d.png", dpi=300)
    plt.close()


def plot_energy_deposition(summary, nominal, cfg, energy, thickness):
    """Overlaid energy-deposition spectra at nominal energy and thickness."""
    values = np.concatenate([nominal[m["id"]]["edep_DUT_MeV"].to_numpy() for m in cfg["materials"]])
    lo, hi = 0.0, float(np.percentile(values, 99.5))
    x_grid = np.linspace(lo, hi, 400)

    fig, ax = plt.subplots(figsize=(8.5, 5.0), dpi=300)
    for mat in cfg["materials"]:
        row = summary[(summary["mat_id"] == mat["id"])
                      & (summary["thickness_mm"] == thickness)
                      & (summary["energy_GeV"] == energy)].iloc[0]
        vals = nominal[mat["id"]]["edep_DUT_MeV"].to_numpy()
        ax.hist(vals, bins=75, range=(lo, hi), density=True, histtype="stepfilled",
                alpha=0.22, color=mat["color"])
        ax.hist(vals, bins=75, range=(lo, hi), density=True, histtype="step",
                linewidth=1.5, color=mat["color"],
                label=f"{mat['label']}  |  MPV = {row['mpv_MeV']:.2f} MeV")
        ax.plot(x_grid, moyal.pdf(x_grid, row["mpv_MeV"], row["moyal_scale"]),
                "--", color=mat["color"], linewidth=1.4)

    ax.set_xlabel(f"Energy Deposited in {thickness:.0f} mm DUT Disc, ΔE [MeV]")
    ax.set_ylabel("Normalized Probability Density [MeV⁻¹]")
    ax.set_title(f"DESY-II {energy:.0f} GeV e⁻ Energy Deposition (S₁ ∩ S₂ Triggered) "
                 f"with Moyal/Landau Fits", weight="bold")
    ax.legend(frameon=True, fontsize=8.5)
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "fig2_landau_edep.png", dpi=300)
    plt.close()


def plot_sipm_yield(summary, optics, mie_eps, cfg, energy, thickness):
    """SiPM photoelectron yield vs. thickness, and Rayleigh vs. Mie collection."""
    opt_cfg = cfg["optical_mc"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0, 4.8), dpi=300,
                                   gridspec_kw={"width_ratios": [1.7, 1.0]})

    for mat in cfg["materials"]:
        sub = summary[(summary["mat_id"] == mat["id"])
                      & (summary["energy_GeV"] == energy)].sort_values("thickness_mm")
        ax1.plot(sub["thickness_mm"], sub["n_pe_nom"], "o-", color=mat["color"],
                 linewidth=2.0, label=mat["label"])
        ax1.fill_between(sub["thickness_mm"], sub["n_pe_min"], sub["n_pe_max"],
                         color=mat["color"], alpha=0.16)
    ax1.axhline(opt_cfg["trigger_threshold_pe"], color="#DC2626", linestyle="--", linewidth=1.5,
                label=f"Clean Trigger Floor ({opt_cfg['trigger_threshold_pe']:.0f} p.e.)")
    ax1.set_xticks(cfg["geometry"]["thickness_scan_mm"])
    ax1.set_xlabel("DUT Disc Thickness [mm]")
    ax1.set_ylabel("Expected SiPM Signal [Photoelectrons / e⁻]")
    ax1.set_title("Hamamatsu S13360-3050CS Signal Yield (PDE = 40%)", weight="bold")
    ax1.legend(frameon=True, fontsize=8)

    bars = [
        ("Mie\n(Method A)", mie_eps, "#94A3B8"),
        ("Rayleigh\n25 wt%", optics[("BiNanoPDMS_25wt", thickness)], "#EA580C"),
        ("Rayleigh\n66 wt%", optics[("BiNanoPDMS_66wt", thickness)], "#DC2626"),
        ("BiPh3-PVK\n(Method C)", optics[("BiPh3PVK", thickness)], "#059669"),
        ("EJ-200\n(Control)", optics[("EJ200", thickness)], "#2563EB"),
    ]
    rects = ax2.bar([b[0] for b in bars], [100.0 * b[1] for b in bars],
                    color=[b[2] for b in bars], edgecolor="k", linewidth=0.6)
    for rect, (_, eps, _) in zip(rects, bars):
        ax2.text(rect.get_x() + rect.get_width() / 2, 100.0 * eps + 0.5, f"{100.0 * eps:.1f}%",
                 ha="center", va="bottom", fontsize=8.5, weight="bold")
    ax2.set_ylabel("Optical Collection Efficiency ε_geo [%]")
    ax2.set_title(f"Optical Collection $\\epsilon_{{\\text{{geo}}}}$\n"
                  f"({thickness:.0f} mm Disc → 3×3 mm² SiPM)", weight="bold", fontsize=11)
    ax2.tick_params(axis="x", labelsize=8, labelrotation=25)

    fig.tight_layout()
    plt.savefig(RESULTS_DIR / "fig3_sipm_yield.png", dpi=300)
    plt.close()


def plot_scattering(summary, cfg, thickness):
    """MCS angle vs. the Highland formula, and Bremsstrahlung rate vs. beam energy."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0, 4.6), dpi=300)
    for mat in cfg["materials"]:
        sub = summary[(summary["mat_id"] == mat["id"])
                      & (summary["thickness_mm"] == thickness)].sort_values("energy_GeV")
        ax1.plot(sub["energy_GeV"], sub["rms_scatter_mrad"], "o-", color=mat["color"],
                 linewidth=2.0, label=f"{mat['label']} (Geant4)")
        ax1.plot(sub["energy_GeV"], sub["highland_mrad"], "--", color=mat["color"],
                 alpha=0.65, linewidth=1.2)
        ax2.plot(sub["energy_GeV"], sub["brem_rate_pct"], "s-", color=mat["color"],
                 linewidth=2.0, label=mat["label"])

    ax1.plot([], [], "k--", label="Analytical Highland Theory")
    ax1.set_xlabel("DESY-II Electron Beam Energy [GeV]")
    ax1.set_ylabel("Planar RMS Scattering Angle θ_RMS [mrad]")
    ax1.set_title(f"Multiple Coulomb Scattering ({thickness:.0f} mm Disc)", weight="bold")
    ax1.legend(frameon=True, fontsize=7.5)

    ax2.set_xlabel("DESY-II Electron Beam Energy [GeV]")
    ax2.set_ylabel("Events with Bremsstrahlung γ > 1 MeV [%]")
    ax2.set_title(f"High-Z Bremsstrahlung Signature ({thickness:.0f} mm Disc)", weight="bold")
    ax2.legend(frameon=True, fontsize=8)

    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "fig4_scattering_brem.png", dpi=300)
    plt.close()


def write_summary_table(summary, optics, mie_eps, mie_ctrl, cfg, energy, thickness):
    """Proposal-ready Markdown benchmark table for the nominal configuration."""
    nominal = summary[(summary["thickness_mm"] == thickness) & (summary["energy_GeV"] == energy)]
    lines = [
        f"# BL4S 2027 Simulation Summary ({cfg['project']['facility']}, "
        f"{energy:.0f} GeV e⁻, {thickness:.0f} mm Disc)\n",
        "| Material | Density (g/cm³) | $X_0$ (mm) | $x/X_0$ (%) | Landau MPV ΔE (MeV) "
        "| Optical $\\epsilon_{\\text{geo}}$ (%) | Expected SiPM Yield (p.e.) "
        "| Geant4 $\\theta_{\\text{RMS}}$ (mrad) | Highland $\\theta_0$ (mrad) | Brem γ >1 MeV (%) |",
        "| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |",
    ]
    for _, r in nominal.iterrows():
        lines.append(
            f"| {r['label']} | {r['density_g_cm3']:.2f} | {r['radlen_mm']:.1f} | {r['x_over_X0_pct']:.2f}% "
            f"| {r['mpv_MeV']:.3f} ± {r['mpv_err_MeV']:.3f} | {r['eps_geo_pct']:.1f}% "
            f"| **{r['n_pe_nom']:.0f}** ({r['n_pe_min']:.0f}–{r['n_pe_max']:.0f}) "
            f"| {r['rms_scatter_mrad']:.2f} | {r['highland_mrad']:.2f} | {r['brem_rate_pct']:.2f}% |"
        )

    rayleigh = optics[("BiNanoPDMS_25wt", thickness)]
    lines.append(
        f"\n*Optical Rayleigh vs. Mie Benchmark ({thickness:.0f} mm disc):* Method A hazy "
        f"microcrystals (`scat_length = {mie_ctrl['scat_length_mm']} mm`) yield "
        f"$\\epsilon_{{\\text{{geo}}}} = {100.0 * mie_eps:.1f}\\%$ vs. Method B transparent Rayleigh "
        f"nanocrystals ($\\epsilon_{{\\text{{geo}}}} = {100.0 * rayleigh:.1f}\\%$), "
        f"a **{rayleigh / max(mie_eps, 1e-6):.1f}×** collection gain."
    )
    (RESULTS_DIR / "summary_table.md").write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--smoke-test", action="store_true", help="fast 2,000-event verification run")
    parser.add_argument("--skip-geant4", action="store_true", help="reuse cached CSVs, re-run optics and plots")
    args = parser.parse_args()

    cfg = json.loads((ROOT / "config.json").read_text())
    RESULTS_DIR.mkdir(exist_ok=True)
    DATA_DIR.mkdir(exist_ok=True)

    events = cfg["project"]["events_smoke" if args.smoke_test else "events_full"]
    nom_E = cfg["beam"]["nominal_energy_GeV"]
    nom_T = cfg["geometry"]["nominal_thickness_mm"]
    jobs = build_jobs(cfg, events)

    if args.skip_geant4:
        print(f"[2/4] Reusing {len(jobs)} cached Geant4 CSVs (--skip-geant4)")
    else:
        build()
        print(f"[2/4] Running {len(jobs)} Geant4 sweeps ({events:,} events each) "
              f"across {cfg['project']['max_workers']} workers...")
        with ProcessPoolExecutor(max_workers=cfg["project"]["max_workers"]) as pool:
            for tag, path in pool.map(run, jobs):
                print(f"   ✓ {tag} -> {path.name}")

    print("[3/4] Running the 3D optical photon Monte Carlo...")
    opt_cfg = cfg["optical_mc"]
    n_photons = 25_000 if args.smoke_test else opt_cfg["n_photons"]
    optics = {
        (mat["id"], t): optical_efficiency(
            cfg, t, mat["rindex"], mat["abs_length_mm"], mat["scat_length_mm"],
            n_photons, cfg["project"]["seed"])
        for mat in cfg["materials"]
        for t in cfg["geometry"]["thickness_scan_mm"]
    }
    mie_ctrl = opt_cfg["mie_control"]
    mie_eps = optical_efficiency(
        cfg, nom_T, mie_ctrl["rindex"], mie_ctrl["abs_length_mm"], mie_ctrl["scat_length_mm"],
        n_photons, cfg["project"]["seed"] + 99)

    print("[4/4] Generating 300-DPI figures and summary_table.md...")
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.25})
    summary, nominal = collect_metrics(cfg, jobs, optics, nom_E, nom_T)
    summary.to_csv(RESULTS_DIR / "all_metrics.csv", index=False)

    plot_setup_3d(cfg)
    plot_energy_deposition(summary, nominal, cfg, nom_E, nom_T)
    plot_sipm_yield(summary, optics, mie_eps, cfg, nom_E, nom_T)
    plot_scattering(summary, cfg, nom_T)
    write_summary_table(summary, optics, mie_eps, mie_ctrl, cfg, nom_E, nom_T)

    print("\nDone! Inspect results/summary_table.md and results/fig1..fig4.png.")


if __name__ == "__main__":
    main()

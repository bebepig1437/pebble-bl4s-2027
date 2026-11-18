"""
Vectorized 3D optical photon Monte Carlo ray-tracer for cylindrical scintillator discs.

Geometry: cylinder (diameter D, thickness L) wrapped in 3M ESR specular foil
(R = 98.5%) except for a centred square SiPM coupled through Dow Corning Q2-3067
grease (n = 1.465). Physics: exact unpolarized Fresnel transmission with total
internal reflection, bulk self-absorption, and isotropic scattering used to
quantify Method B transparent Rayleigh nanocrystals against Method A hazy Mie
microcrystals.

Returns the optical collection efficiency eps_geo = (photons reaching the SiPM
aperture and surviving both interfaces) / (photons generated).
"""

import numpy as np


def _fresnel_transmission(cos_i, n1, n2):
    """Unpolarized Fresnel transmission coefficient from medium n1 to n2 (handles TIR)."""
    cos_i = np.clip(cos_i, 0.0, 1.0)
    sin_t2 = (n1 / n2) ** 2 * (1.0 - cos_i * cos_i)

    trans = np.zeros_like(cos_i)
    transmitted = sin_t2 < 1.0
    if np.any(transmitted):
        ci = cos_i[transmitted]
        ct = np.sqrt(1.0 - sin_t2[transmitted])
        rs = ((n1 * ci - n2 * ct) / (n1 * ci + n2 * ct)) ** 2
        rp = ((n1 * ct - n2 * ci) / (n1 * ct + n2 * ci)) ** 2
        trans[transmitted] = 1.0 - 0.5 * (rs + rp)
    return trans


def simulate_optical_collection(
    disc_diameter_mm,
    thickness_mm,
    trigger_xy_mm,
    sipm_xy_mm,
    rindex,
    abs_length_mm,
    scat_length_mm,
    grease_rindex=1.465,
    sipm_epoxy_rindex=1.55,
    esr_reflectivity=0.985,
    n_photons=100_000,
    max_bounces=200,
    seed=2027,
):
    """Ray-trace isotropic photons and return the optical collection efficiency."""
    rng = np.random.default_rng(seed)
    R = 0.5 * disc_diameter_mm
    half_L = 0.5 * thickness_mm
    half_trig = 0.5 * trigger_xy_mm
    half_sipm = 0.5 * sipm_xy_mm

    # Photons are generated uniformly in the S1^S2 trigger region of the disc.
    pos = np.empty((n_photons, 3))
    pos[:, 0] = rng.uniform(-half_trig, half_trig, n_photons)
    pos[:, 1] = rng.uniform(-half_trig, half_trig, n_photons)
    pos[:, 2] = rng.uniform(-half_L, half_L, n_photons)

    cos_t = rng.uniform(-1.0, 1.0, n_photons)
    sin_t = np.sqrt(np.maximum(0.0, 1.0 - cos_t * cos_t))
    phi = rng.uniform(0.0, 2.0 * np.pi, n_photons)
    dirs = np.column_stack((sin_t * np.cos(phi), sin_t * np.sin(phi), cos_t))

    detected = 0
    eps = 1e-6

    for _ in range(max_bounces):
        if len(pos) == 0:
            break

        x, y, z = pos[:, 0], pos[:, 1], pos[:, 2]
        dx, dy, dz = dirs[:, 0], dirs[:, 1], dirs[:, 2]

        # Ray intersections with the side wall (x^2 + y^2 = R^2) and both caps.
        # Division by zero yields inf/nan, which the masks below discard.
        a = dx * dx + dy * dy
        b = 2.0 * (x * dx + y * dy)
        c = x * x + y * y - R * R
        with np.errstate(invalid="ignore", divide="ignore"):
            t_side = (-b + np.sqrt(np.maximum(0.0, b * b - 4.0 * a * c))) / (2.0 * a)
            t_back = (half_L - z) / dz
            t_front = (-half_L - z) / dz
        t_side = np.where((a > 1e-12) & (t_side > 0.0), t_side, np.inf)
        t_back = np.where(dz > 1e-12, t_back, np.inf)
        t_front = np.where(dz < -1e-12, t_front, np.inf)
        t_wall = np.minimum(t_side, np.minimum(t_back, t_front))

        # Competing bulk absorption and scattering free paths.
        s_abs = rng.exponential(abs_length_mm, len(pos))
        s_scat = rng.exponential(scat_length_mm, len(pos))
        is_abs = (s_abs < t_wall) & (s_abs <= s_scat)
        is_scat = (s_scat < t_wall) & (s_scat < s_abs)
        is_wall = ~(is_abs | is_scat)

        next_pos, next_dir = [], []

        # Bulk scattering: advance to the scatter point, then re-isotropize.
        if np.any(is_scat):
            p_sc = pos[is_scat] + dirs[is_scat] * s_scat[is_scat, None]
            ct = rng.uniform(-1.0, 1.0, len(p_sc))
            st = np.sqrt(np.maximum(0.0, 1.0 - ct * ct))
            ph = rng.uniform(0.0, 2.0 * np.pi, len(p_sc))
            next_pos.append(p_sc)
            next_dir.append(np.column_stack((st * np.cos(ph), st * np.sin(ph), ct)))

        if np.any(is_wall):
            pw = pos[is_wall] + dirs[is_wall] * np.maximum(0.0, t_wall[is_wall] - eps)[:, None]
            dw = dirs[is_wall].copy()

            hit_back = t_back[is_wall] <= np.minimum(t_front[is_wall], t_side[is_wall]) + 1e-9
            hit_front = (~hit_back) & (t_front[is_wall] <= t_side[is_wall] + 1e-9)
            hit_side = ~(hit_back | hit_front)

            # SiPM aperture on the back face: Fresnel through DUT -> grease -> epoxy.
            on_sipm = hit_back & (np.abs(pw[:, 0]) <= half_sipm) & (np.abs(pw[:, 1]) <= half_sipm)
            survive = np.zeros(len(pw), dtype=bool)

            if np.any(on_sipm):
                cos_i = np.abs(dw[on_sipm, 2])
                t1 = _fresnel_transmission(cos_i, rindex, grease_rindex)
                sin_g2 = np.minimum(1.0, (rindex / grease_rindex) ** 2 * (1.0 - cos_i * cos_i))
                cos_g = np.sqrt(np.maximum(0.0, 1.0 - sin_g2))
                t2 = _fresnel_transmission(cos_g, grease_rindex, sipm_epoxy_rindex)

                sipm_idx = np.where(on_sipm)[0]
                detected_now = rng.random(len(sipm_idx)) < t1 * t2
                detected += int(np.count_nonzero(detected_now))

                # Photons reflected off the SiPM window bounce back into the disc.
                reflected = sipm_idx[~detected_now]
                dw[reflected, 2] *= -1.0
                survive[reflected] = True

            # Every non-SiPM face is wrapped in 3M ESR specular foil.
            on_esr = ~on_sipm
            if np.any(on_esr):
                esr_idx = np.where(on_esr)[0]
                kept = esr_idx[rng.random(len(esr_idx)) < esr_reflectivity]
                survive[kept] = True

                flat = kept[hit_front[kept] | hit_back[kept]]
                dw[flat, 2] *= -1.0

                side = kept[hit_side[kept]]
                if len(side) > 0:
                    nx = pw[side, 0] / R
                    ny = pw[side, 1] / R
                    dot = dw[side, 0] * nx + dw[side, 1] * ny
                    dw[side, 0] -= 2.0 * dot * nx
                    dw[side, 1] -= 2.0 * dot * ny

            if np.any(survive):
                next_pos.append(pw[survive])
                next_dir.append(dw[survive])

        if not next_pos:
            break
        pos = np.vstack(next_pos)
        dirs = np.vstack(next_dir)

    return detected / n_photons

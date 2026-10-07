import sys
import numpy as np

def load_profile(path, Lz, dx=0.1):
    data = np.load(path, allow_pickle=True)

    if data.dtype == object:
        data = data.reshape(-1)[0]

        if isinstance(data, dict):
            return data["z"], data["com"]

    edges = np.linspace(-Lz / 2, Lz / 2, int(round(Lz / dx)) + 1)
    counts, _ = np.histogram(data, bins=edges)

    z = 0.5 * (edges[1:] + edges[:-1])

    return z, counts.astype(float)

def build_bias_table(profile_path, table_path, T, Lz, old_table_path=None, dx=0.1):
    z, c = load_profile(profile_path, Lz, dx)

    xg = np.linspace(-Lz / 2, Lz / 2, int(round(Lz / dx)) + 1) # create grid
    cg = np.interp(xg, z, c)

    # Smooth the profile.
    sigma = 0.3
    k = np.arange(-int(4 * sigma / dx), int(4 * sigma / dx) + 1)
    w = np.exp(-0.5 * (k * dx / sigma) ** 2)
    w /= w.sum()

    cg = np.convolve(np.pad(cg, len(k) // 2, mode="edge"), w, mode="valid")

    # Bilayer symmetry
    cg = 0.5 * (cg + cg[::-1])

    # Avoid log(0)
    cfloor = 1e-2 * cg.max()
    sampled = cg > cfloor
    cg = np.maximum(cg, cfloor)

    RT = 0.008314463 * T

    if old_table_path is None:
        u = RT * np.log(cg)
    else:
        u_old = np.load(old_table_path)

        c_target = cg[sampled].mean()
        u = u_old + RT * np.log(cg / c_target)

    # Set zero in bulk water and enforce periodicity.
    u -= u[0]
    u = np.clip(u, -25.0, 25.0)
    u[-1] = u[0]

    np.save(table_path, u)

profile_path = sys.argv[1]
table_path = sys.argv[2]
T = 323.15
Lz = 22.0

old_table_path = None
if len(sys.argv) == 4:
    old_table_path = sys.argv[3]

build_bias_table(profile_path, table_path, T, Lz, old_table_path=old_table_path)
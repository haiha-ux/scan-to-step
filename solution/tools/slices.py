import numpy as np, sys, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, 'solution')
import meas as M
i = int(sys.argv[1])
m, mm, to, good = M.load(i)
V = mm.vertices[np.unique(mm.faces[good].ravel())]
specs = [x.split(':') for x in sys.argv[2:]]   # axis:lo:hi:u:v
fig, axs = plt.subplots(1, len(specs), figsize=(7 * len(specs), 7))
axs = np.atleast_1d(axs)
for ax_, (a, lo, hi, u, v) in zip(axs, specs):
    a, u, v = int(a), int(u), int(v); lo, hi = float(lo), float(hi)
    s = (V[:, a] > lo) & (V[:, a] < hi)
    ax_.scatter(V[s, u], V[s, v], s=0.2)
    ax_.set_aspect('equal'); ax_.grid(True, lw=0.3); ax_.set_title(f'axis{a} in [{lo},{hi}]  ({"xyz"[u]},{"xyz"[v]})')
fig.tight_layout(); fig.savefig(f'solution/img/_slice_{i}.png', dpi=80)

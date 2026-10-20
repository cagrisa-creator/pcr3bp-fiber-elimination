#!/usr/bin/env python3
"""
Regenerates paper/figures/phase_map.png from results/phase_map_data.csv.

v1.3.0: the horizontal axis is m_r, the mass fraction of the REGULARIZED
primary (m_r = 1 - mu, with mu the distant-primary mass parameter of the
defining function). Kim's 16/17 threshold is drawn as a reference line: his
theorem guarantees non-convexity of the heavy-primary component (m_r > 1/2)
for m_r < 16/17 near the critical energy; the reconnaissance shows the
non-convex window in fact extends well beyond it on the heavy side.

Author: Cesar Augusto Grisa (UFF, Brazil).
"""
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

rows = list(csv.DictReader(open('../results/phase_map_data.csv')))
mr  = [float(r['m_r']) for r in rows]
dc  = [float(r['delta_c']) for r in rows]
D   = [float(r['D_true_min']) for r in rows]

fig, ax = plt.subplots(figsize=(8.6, 4.6))
for m, d, v in zip(mr, dc, D):
    ax.scatter(m, d, c=('tab:red' if v < 0 else 'tab:green'),
               s=42, edgecolors='k', linewidths=0.4, zorder=3)
ax.set_yscale('log')
ax.set_xlabel(r'$m_r$ (mass fraction of the regularized primary)')
ax.set_ylabel(r'$\Delta c = c - c_{L_1}$')
ax.axvline(16/17, color='0.35', lw=1.0, ls='--', zorder=2)
ax.text(16/17, 1.6, r"Kim $16/17$", ha='center', fontsize=8, color='0.25')
ax.axvline(0.987849414390376, color='tab:blue', lw=0.9, ls=':', zorder=2)
ax.text(0.9878, 2.6, "Earth", ha='center', fontsize=8, color='tab:blue')
ax.axvline(0.012150585609624, color='tab:blue', lw=0.9, ls=':', zorder=2)
ax.text(0.0122, 2.6, "Moon", ha='center', fontsize=8, color='tab:blue')
ax.set_xlim(-0.03, 1.03)
ax.set_ylim(5e-5, 4.5)
ax.set_title('LC-convexity reconnaissance: green = both gates positive, '
             'red = determinant gate negative', fontsize=10)
fig.tight_layout()
fig.savefig('../paper/figures/phase_map.png', dpi=150)
print('wrote ../paper/figures/phase_map.png')

#!/usr/bin/env python3
"""
Virtual Screening Binding Energy Distribution
Run on WSL2:
    cd ~/projects/rtc-pan-coronavirus
    python3 CGCP/scripts/phase-2/screening_curve.py
"""

import os, glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MultipleLocator

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'figure.facecolor': 'white',
    'axes.facecolor':   'white',
    'axes.spines.top':  False,
    'axes.spines.right':False,
    'axes.grid':        False,
    'xtick.direction':  'out',
    'ytick.direction':  'out',
    'xtick.major.width':1.2,
    'ytick.major.width':1.2,
    'xtick.major.size': 5,
    'ytick.major.size': 5,
    'pdf.fonttype': 42,
})

BASE    = os.path.expanduser('~/projects/rtc-pan-coronavirus')
RESULTS = f'{BASE}/04-hits/CGCP_run2_tier2_10M'
OUT     = f'{BASE}/04-hits'

# ── Load all scores ───────────────────────────────────────────
print("Loading scores...", flush=True)
scores = []
for f in sorted(glob.glob(f'{RESULTS}/scores_task*.tsv')):
    with open(f) as fh:
        for line in fh:
            if line.startswith('compound') or line.startswith('ERROR'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                try:
                    scores.append(float(parts[2]))
                except:
                    pass

scores = np.array(scores)
print(f"  {len(scores):,} compounds loaded")
print(f"  Range: {scores.min():.3f} to {scores.max():.3f} kcal/mol")
print(f"  Mean:  {scores.mean():.3f}  Median: {np.median(scores):.3f}")

# ── Thresholds ────────────────────────────────────────────────
T1, T2, T3 = -7.0, -8.0, -9.0
n_t1 = (scores <= T1).sum()
n_t2 = (scores <= T2).sum()
n_t3 = (scores <= T3).sum()

print(f"\n  ≤ {T1}: {n_t1:,}  ({n_t1/len(scores)*100:.2f}%)")
print(f"  ≤ {T2}: {n_t2:,}  ({n_t2/len(scores)*100:.2f}%)")
print(f"  ≤ {T3}: {n_t3:,}  ({n_t3/len(scores)*100:.3f}%)")

# ════════════════════════════════════════════════════════════
# FIGURE — 2 panels: histogram + cumulative curve
# ════════════════════════════════════════════════════════════
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))
fig.patch.set_facecolor('white')
fig.suptitle(
    'Virtual Screening — Binding Energy Distribution\n'
    'NSP12–NSP7 Interface  |  Enamine REAL Tier2  |  AutoDock Vina',
    fontsize=13, fontweight='bold', y=1.02)

# ── Colour zones ──────────────────────────────────────────────
C_WEAK   = '#CCCCCC'
C_MOD    = '#7FBADC'
C_STRONG = '#1A5276'
C_TOP    = '#C0392B'

bins = np.arange(np.floor(scores.min()) - 0.5, scores.max() + 1.0, 0.25)

# ── Panel 1: Histogram ────────────────────────────────────────
ax1.set_facecolor('white')
counts, edges, patches = ax1.hist(scores, bins=bins,
                                   color=C_WEAK, edgecolor='white',
                                   linewidth=0.3, zorder=3)

# Colour zones
for patch, left in zip(patches, edges[:-1]):
    if left <= T3:
        patch.set_facecolor(C_TOP)
    elif left <= T2:
        patch.set_facecolor(C_STRONG)
    elif left <= T1:
        patch.set_facecolor(C_MOD)

# Threshold lines
for thresh, col, lbl in [
        (T1, C_MOD,    f'≤ {T1} ({n_t1:,})'),
        (T2, C_STRONG, f'≤ {T2} ({n_t2:,})'),
        (T3, C_TOP,    f'≤ {T3} ({n_t3})')]:
    ax1.axvline(thresh, color=col, lw=1.8, ls='--', zorder=5, alpha=0.9)
    ax1.text(thresh - 0.05, ax1.get_ylim()[1] * 0.97 if ax1.get_ylim()[1] > 0 else 1,
             lbl, ha='right', va='top', fontsize=8,
             color=col, fontweight='bold', rotation=90)

ax1.set_xlabel('Binding energy (kcal/mol)', fontsize=11, fontweight='bold')
ax1.set_ylabel('Number of compounds', fontsize=11, fontweight='bold')
ax1.set_title('Binding Energy Histogram', fontsize=11, fontweight='bold', pad=8)
ax1.xaxis.set_minor_locator(MultipleLocator(0.5))
ax1.spines['left'].set_position(('outward', 5))
ax1.spines['bottom'].set_position(('outward', 5))

# Stats box
stats_txt = (f"n = {len(scores):,}\n"
             f"Mean = {scores.mean():.2f}\n"
             f"Median = {np.median(scores):.2f}\n"
             f"Best = {scores.min():.3f}")
ax1.text(0.97, 0.97, stats_txt, transform=ax1.transAxes,
         ha='right', va='top', fontsize=8.5, color='#333333',
         bbox=dict(boxstyle='round,pad=0.4', facecolor='#F4F8FB',
                   edgecolor='#CCCCCC', linewidth=0.8))

# ── Panel 2: Cumulative distribution ─────────────────────────
sorted_scores = np.sort(scores)
cum_pct       = np.arange(1, len(sorted_scores)+1) / len(sorted_scores) * 100

ax2.set_facecolor('white')

# Shade regions
for (xlo, xhi, col, alpha) in [
        (sorted_scores.min(), T1, '#EAF4FB', 0.6),
        (T1, T2, '#D6EAF8', 0.6),
        (T2, T3, '#1A5276', 0.12),
        (T3, sorted_scores.max(), C_TOP, 0.18)]:
    ax2.axvspan(xlo, xhi, color=col, alpha=alpha, zorder=1)

ax2.plot(sorted_scores, cum_pct,
         color='#0D2137', lw=2.0, zorder=4)

# Threshold markers
for thresh, col, n in [
        (T1, C_MOD,    n_t1),
        (T2, C_STRONG, n_t2),
        (T3, C_TOP,    n_t3)]:
    pct = n / len(scores) * 100
    ax2.axvline(thresh, color=col, lw=1.8, ls='--', zorder=5)
    ax2.plot(thresh, pct, 'o', color=col, ms=7, zorder=6)
    ax2.annotate(f'{pct:.2f}%\n({n:,})',
                 xy=(thresh, pct),
                 xytext=(thresh + 0.15, pct + 3),
                 fontsize=8, color=col, fontweight='bold',
                 arrowprops=dict(arrowstyle='->', color=col, lw=1.0))

ax2.set_xlabel('Binding energy (kcal/mol)', fontsize=11, fontweight='bold')
ax2.set_ylabel('Cumulative compounds (%)', fontsize=11, fontweight='bold')
ax2.set_title('Cumulative Distribution', fontsize=11, fontweight='bold', pad=8)
ax2.set_ylim(0, 102)
ax2.spines['left'].set_position(('outward', 5))
ax2.spines['bottom'].set_position(('outward', 5))

# Legend
legend_items = [
    mpatches.Patch(color=C_WEAK,   label=f'> {T1} kcal/mol'),
    mpatches.Patch(color=C_MOD,    label=f'{T2} to {T1} kcal/mol  ({n_t1-n_t2:,})'),
    mpatches.Patch(color=C_STRONG, label=f'{T3} to {T2} kcal/mol  ({n_t2-n_t3:,})'),
    mpatches.Patch(color=C_TOP,    label=f'≤ {T3} kcal/mol  ({n_t3})  ← top hits'),
]
fig.legend(handles=legend_items, loc='lower center', ncol=4,
           fontsize=9, frameon=True, facecolor='white',
           edgecolor='#CCCCCC', bbox_to_anchor=(0.5, -0.06))

plt.tight_layout(rect=[0, 0.06, 1, 1])

for ext in ('png', 'pdf'):
    p = f'{OUT}/Fig_ScreeningDistribution.{ext}'
    fig.savefig(p, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved: {p}")

plt.close()

#!/usr/bin/env python3
"""
Script 16 - Hit Analysis + 3D Visualization
Pan-coronavirus RTC Inhibitor Discovery Pipeline
Olivier Nsekuye | GIGA-VIN Lab, University of Liege | 2026-03-10

Identifies hits, dual-target compounds, generates score distributions,
and prepares 3D visualization of top hits in binding pockets.
"""

import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
from datetime import datetime

BASE    = Path(__file__).parent.parent
SCORES  = BASE / "04-hits" / "scores"
PDBQT   = BASE / "data" / "pdbqt"
FIGS    = BASE / "figures" / "script16"
HITS    = BASE / "04-hits"
FIGS.mkdir(parents=True, exist_ok=True)

TARGETS = ["NSP12-NSP7", "NSP9-NSP12", "NSP12-NSP8"]
HIT_THRESHOLD   = -7.0   # kcal/mol
STRONG_THRESHOLD = -8.0  # kcal/mol

COLORS = {
    "NSP12-NSP7": "#2196F3",
    "NSP9-NSP12": "#F44336",
    "NSP12-NSP8": "#4CAF50",
}

# ── Load scores ───────────────────────────────────────────────────────────
def load_scores():
    dfs = {}
    for target in TARGETS:
        f = SCORES / f"all_scores_{target}.tsv"
        df = pd.read_csv(f, sep='\t', header=None,
                         names=['zinc_id', 'score'])
        df = df[df['score'] != 'ERROR']
        df['score'] = pd.to_numeric(df['score'], errors='coerce')
        df = df.dropna(subset=['score'])
        df['zinc_id'] = df['zinc_id'].astype(str)
        df = df.sort_values('score')
        dfs[target] = df
        print(f"  {target}: {len(df):,} compounds loaded")
    return dfs

# ── Hit identification ────────────────────────────────────────────────────
def identify_hits(dfs):
    hits = {}
    for target, df in dfs.items():
        h = df[df['score'] <= HIT_THRESHOLD].copy()
        h['target'] = target
        hits[target] = h
        strong = len(h[h['score'] <= STRONG_THRESHOLD])
        print(f"  {target}: {len(h):,} hits (<= {HIT_THRESHOLD}) "
              f"| {strong} strong (<= {STRONG_THRESHOLD})")
    return hits

# ── Dual/triple target hits ───────────────────────────────────────────────
def find_multitarget_hits(hits):
    hit_sets = {t: set(h['zinc_id']) for t, h in hits.items()}

    dual_12_7_9_12  = hit_sets["NSP12-NSP7"] & hit_sets["NSP9-NSP12"]
    dual_12_7_12_8  = hit_sets["NSP12-NSP7"] & hit_sets["NSP12-NSP8"]
    dual_9_12_12_8  = hit_sets["NSP9-NSP12"] & hit_sets["NSP12-NSP8"]
    triple          = hit_sets["NSP12-NSP7"] & hit_sets["NSP9-NSP12"] & hit_sets["NSP12-NSP8"]

    print(f"\n  Dual-target hits:")
    print(f"    NSP12-NSP7 + NSP9-NSP12:  {len(dual_12_7_9_12)}")
    print(f"    NSP12-NSP7 + NSP12-NSP8:  {len(dual_12_7_12_8)}")
    print(f"    NSP9-NSP12 + NSP12-NSP8:  {len(dual_9_12_12_8)}")
    print(f"  Triple-target hits:          {len(triple)}")

    return {
        "dual_NSP12-NSP7_NSP9-NSP12": sorted(dual_12_7_9_12),
        "dual_NSP12-NSP7_NSP12-NSP8": sorted(dual_12_7_12_8),
        "dual_NSP9-NSP12_NSP12-NSP8": sorted(dual_9_12_12_8),
        "triple":                      sorted(triple),
    }

# ── Score distribution plot ───────────────────────────────────────────────
def plot_score_distributions(dfs):
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle("AutoDock Vina Score Distributions — Pan-CoV RTC Screen\n"
                 "ZINC20 Drug-Like Subset (9,800 compounds per target)",
                 fontsize=13, fontweight='bold', y=1.02)

    for ax, target in zip(axes, TARGETS):
        df    = dfs[target]
        color = COLORS[target]
        scores = df['score'].values
        n_hits   = (scores <= HIT_THRESHOLD).sum()
        n_strong = (scores <= STRONG_THRESHOLD).sum()

        ax.hist(scores, bins=60, color=color, alpha=0.75, edgecolor='white', lw=0.3)
        ax.axvline(HIT_THRESHOLD,   color='orange', lw=2, ls='--',
                   label=f'Hit threshold ({HIT_THRESHOLD})')
        ax.axvline(STRONG_THRESHOLD, color='red',    lw=2, ls='--',
                   label=f'Strong hit ({STRONG_THRESHOLD})')

        ax.set_xlabel('Vina Score (kcal/mol)', fontsize=11)
        ax.set_ylabel('Count', fontsize=11)
        ax.set_title(f'{target}\nn={len(df):,} | hits={n_hits} | strong={n_strong}',
                     fontsize=11, fontweight='bold')
        ax.legend(fontsize=8)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.tight_layout()
    out = FIGS / "score_distributions.png"
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out}")

# ── Top hits table ────────────────────────────────────────────────────────
def build_top_hits_table(dfs, hits, multitarget):
    triple_set = set(multitarget["triple"])
    dual_sets  = {
        "NSP12-NSP7+NSP9-NSP12": set(multitarget["dual_NSP12-NSP7_NSP9-NSP12"]),
        "NSP12-NSP7+NSP12-NSP8": set(multitarget["dual_NSP12-NSP7_NSP12-NSP8"]),
        "NSP9-NSP12+NSP12-NSP8": set(multitarget["dual_NSP9-NSP12_NSP12-NSP8"]),
    }

    rows = []
    for target, df in dfs.items():
        top = df.head(20)
        for _, row in top.iterrows():
            zid = str(row['zinc_id'])
            score = row['score']
            if zid in triple_set:
                activity = "TRIPLE"
            else:
                duals = [k for k, s in dual_sets.items() if zid in s]
                activity = duals[0] if duals else "single"
            rows.append({
                "zinc_id":  zid,
                "target":   target,
                "score":    score,
                "activity": activity,
            })

    table = pd.DataFrame(rows).sort_values('score')
    out   = HITS / "top_hits_table.tsv"
    table.to_csv(out, sep='\t', index=False)
    print(f"  Saved: {out}")
    return table

# ── Multitarget scatter plot ──────────────────────────────────────────────
def plot_multitarget_scatter(dfs, multitarget):
    df1 = dfs["NSP12-NSP7"].set_index('zinc_id')
    df2 = dfs["NSP12-NSP8"].set_index('zinc_id')
    common = df1.index.intersection(df2.index)
    merged = pd.DataFrame({
        'NSP12-NSP7': df1.loc[common, 'score'],
        'NSP12-NSP8': df2.loc[common, 'score'],
    })

    dual_set   = set(multitarget["dual_NSP12-NSP7_NSP12-NSP8"])
    triple_set = set(multitarget["triple"])

    fig, ax = plt.subplots(figsize=(8, 7))

    # All compounds
    ax.scatter(merged['NSP12-NSP7'], merged['NSP12-NSP8'],
               alpha=0.15, s=8, color='#BBBBBB', label='All compounds')

    # Dual hits
    dual_data = merged[merged.index.isin(dual_set)]
    ax.scatter(dual_data['NSP12-NSP7'], dual_data['NSP12-NSP8'],
               alpha=0.8, s=40, color='#FF9800', zorder=3,
               label=f'Dual hits (n={len(dual_data)})')

    # Triple hits
    if triple_set:
        triple_data = merged[merged.index.isin(triple_set)]
        ax.scatter(triple_data['NSP12-NSP7'], triple_data['NSP12-NSP8'],
                   alpha=1.0, s=80, color='#E91E63', zorder=4,
                   marker='*', label=f'Triple hits (n={len(triple_data)})')

    ax.axvline(HIT_THRESHOLD, color='blue',  ls='--', lw=1, alpha=0.5)
    ax.axhline(HIT_THRESHOLD, color='green', ls='--', lw=1, alpha=0.5)
    ax.set_xlabel('NSP12-NSP7 Score (kcal/mol)', fontsize=12)
    ax.set_ylabel('NSP12-NSP8 Score (kcal/mol)', fontsize=12)
    ax.set_title('Dual-Target Activity: NSP12-NSP7 vs NSP12-NSP8\n'
                 'Pan-CoV RTC Virtual Screen', fontsize=12, fontweight='bold')
    ax.legend(fontsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    out = FIGS / "dual_target_scatter.png"
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out}")

# ── Main ──────────────────────────────────────────────────────────────────
def main():
    print("="*65)
    print("Script 16 - Hit Analysis + Visualization")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*65)

    print("\n[1] Loading scores...")
    dfs = load_scores()

    print("\n[2] Identifying hits...")
    hits = identify_hits(dfs)

    print("\n[3] Finding multi-target hits...")
    multitarget = find_multitarget_hits(hits)

    print("\n[4] Plotting score distributions...")
    plot_score_distributions(dfs)

    print("\n[5] Building top hits table...")
    table = build_top_hits_table(dfs, hits, multitarget)

    print("\n[6] Plotting dual-target scatter...")
    plot_multitarget_scatter(dfs, multitarget)

    # Save summary JSON
    summary = {
        "completed": datetime.now().isoformat(),
        "hit_threshold": HIT_THRESHOLD,
        "strong_threshold": STRONG_THRESHOLD,
        "per_target": {
            t: {
                "total": len(dfs[t]),
                "hits":  int((dfs[t]['score'] <= HIT_THRESHOLD).sum()),
                "strong": int((dfs[t]['score'] <= STRONG_THRESHOLD).sum()),
                "best_score": float(dfs[t]['score'].min()),
                "best_zinc":  str(dfs[t].iloc[0]['zinc_id']),
            } for t in TARGETS
        },
        "multitarget": {k: len(v) for k, v in multitarget.items()},
        "multitarget_ids": multitarget,
    }
    out = HITS / "hit_analysis_summary.json"
    with open(out, 'w') as f:
        json.dump(summary, f, indent=2)

    print("\n" + "="*65)
    print("HIT ANALYSIS SUMMARY")
    print("="*65)
    for t in TARGETS:
        s = summary["per_target"][t]
        print(f"  {t:<20} best={s['best_score']:.3f}  "
              f"hits={s['hits']}  strong={s['strong']}  "
              f"top={s['best_zinc']}")
    print(f"\n  Multi-target:")
    for k, v in summary["multitarget"].items():
        print(f"    {k:<40} {v}")
    print(f"\nOutput: {HITS}")
    print("Status: COMPLETE")

if __name__ == "__main__":
    main()

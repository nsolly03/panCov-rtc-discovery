import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

categories = ['Anchor', 'Aromatic', 'Hydrophobic', 'Charged+', 'Charged−', 'H-bond']
values     = [1, 1, 40, 8, 5, 20]
colors     = ['#555555', '#C0392B', '#95A5A6', '#5B8DB8', '#82C48A', '#C9A84C']
val_colors = ['#555555', '#C0392B', '#7F8C8D', '#2471A3', '#1E8449', '#9A7D0A']

fig, ax = plt.subplots(figsize=(8, 5))
fig.patch.set_facecolor('white')
ax.set_facecolor('white')

x = np.arange(len(categories))
bars = ax.bar(x, values, width=0.42, color=colors, linewidth=0.8,
              edgecolor=[c for c in colors])

# Value labels in matching colors
for bar, val, vc in zip(bars, values, val_colors):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
            str(val), ha='center', va='bottom', fontsize=12,
            fontweight='bold', color=vc, fontfamily='Arial')

# Prism axes: bottom + left only, ending at last tick
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(1.0)
ax.spines['bottom'].set_linewidth(1.0)

# Outward ticks
ax.tick_params(axis='both', direction='out', length=4, width=1.0)
ax.tick_params(axis='x', bottom=True, top=False)
ax.tick_params(axis='y', left=True, right=False)

# Y axis: 0–45, ticks every 5
ax.set_ylim(0, 45)
ax.yaxis.set_major_locator(ticker.MultipleLocator(5))

# Labels
ax.set_xticks(x)
ax.set_xticklabels(categories, fontfamily='Arial', fontsize=11)
ax.set_ylabel('Number of residues', fontfamily='Arial', fontsize=12, fontweight='bold')
ax.set_xlabel('Residue class', fontfamily='Arial', fontsize=12, fontweight='bold')

# Stop axes at last tick
ax.spines['left'].set_bounds(0, 40)
ax.spines['bottom'].set_bounds(x[0] - 0.3, x[-1] + 0.3)

plt.tight_layout()
plt.savefig('interface_residue_classes.pdf', dpi=300, bbox_inches='tight')
plt.savefig('interface_residue_classes.png', dpi=300, bbox_inches='tight')
plt.show()
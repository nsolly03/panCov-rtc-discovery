"""
CGCP PyMOL Publication-Quality 3D Visualization Template
=========================================================
Generates all 3D figures for one interface following the
style of the JB12/JB16.1 presentation but with:
  - Cluster pharmacophore overlay
  - Conservation coloring
  - Positive + negative control panels
  - Consistent Prism-compatible color scheme

Usage:
    pymol -c pymol_cgcp_visuals.py -- NSP12-NSP7 PHE440

Or interactively:
    run pymol_cgcp_visuals.py
    generate_all_views("NSP12-NSP7", "PHE440", chain_a="A", chain_b="C")
"""

import pymol
from pymol import cmd
import os
import sys

# ─────────────────────────────────────────────────────────────
# PRISM-COMPATIBLE COLOR PALETTE
# Matches matplotlib palette used in Phase 0 plots
# ─────────────────────────────────────────────────────────────
COLORS = {
    # Chain colors
    "chain_a":     [0.40, 0.53, 0.67],   # steel blue  #2166AC
    "chain_b":     [0.30, 0.43, 0.07],   # forest green #4DAC26
    "chain_bg":    [0.85, 0.85, 0.85],   # light gray for background chain

    # Hotspot tiers
    "anchor":      [0.84, 0.10, 0.11],   # crimson      #D7191C — primary anchor
    "cluster":     [0.94, 0.49, 0.15],   # orange       #F07B26 — cluster members
    "conserved":   [0.18, 0.37, 0.04],   # dark green   #2E5E09 — cons >= 0.8
    "variable":    [0.84, 0.10, 0.11],   # red          #D7191C — cons < 0.5
    "moderate":    [0.95, 0.60, 0.00],   # amber        #F29900 — cons 0.5-0.8

    # Ligand / pharmacophore
    "ligand":      [1.00, 0.90, 0.00],   # yellow       — ligand sticks
    "pharmacophore": [0.55, 0.00, 0.78], # purple       — pharmacophore spheres
    "hbond":       [0.20, 0.60, 1.00],   # blue dashed  — hydrogen bonds
    "saltbridge":  [1.00, 0.10, 0.10],   # red dashed   — salt bridges

    # Surface
    "surface_main": [0.73, 0.81, 0.89],  # light blue   — primary chain surface
    "surface_alt":  [0.78, 0.91, 0.78],  # light green  — secondary chain surface
}


def setup_pymol_prism():
    """
    Apply Prism-style rendering settings.
    White background, ray-traced, no fog, clean lighting.
    Call once at start of session.
    """
    cmd.set("bg_rgb", [1, 1, 1])             # white background
    cmd.set("ray_opaque_background", 0)       # transparent option
    cmd.set("orthoscopic", 1)                 # no perspective distortion
    cmd.set("antialias", 2)                   # smooth edges
    cmd.set("ray_shadows", 0)                 # clean, no shadows
    cmd.set("ray_trace_fog", 0)               # no fog
    cmd.set("depth_cue", 0)                   # no depth fading
    cmd.set("hash_max", 300)                  # memory for ray tracing
    cmd.set("ray_trace_mode", 1)              # natural + black outline
    cmd.set("cartoon_fancy_helices", 1)       # Prism-style helices
    cmd.set("cartoon_fancy_sheets", 1)        # Prism-style sheets
    cmd.set("cartoon_tube_radius", 0.20)      # thinner tubes
    cmd.set("cartoon_loop_radius", 0.15)      # thinner loops
    cmd.set("stick_radius", 0.12)             # thin sticks like presentation
    cmd.set("sphere_scale", 0.35)             # compact spheres for hotspots
    cmd.set("label_size", 10)                 # readable labels
    cmd.set("label_color", "black")
    cmd.set("label_font_id", 7)               # Arial equivalent


def render_figure(filename, width=2400, height=1800, dpi=300):
    """
    Ray-trace and save at publication resolution.
    Matches 8x6 inch @ 300 DPI = Nature/Cell standard.
    """
    cmd.ray(width, height, renderer=0)
    cmd.png(filename, dpi=dpi, quiet=1)
    print(f"  Saved: {filename}")


def color_by_conservation(selection, conservation_dict):
    """
    Color residues by conservation score using Prism palette.
    conservation_dict: {residue_number: score_0_to_1}

    Color scheme:
        1.000      = dark green  (pan-coronavirus)
        0.8-0.999  = green
        0.5-0.799  = amber
        < 0.5      = red (SARS-selective or variable)
    """
    for resi, score in conservation_dict.items():
        sel = f"({selection}) and resi {resi}"
        if score >= 0.999:
            cmd.color("forest",   sel)
        elif score >= 0.800:
            cmd.color("green",    sel)
        elif score >= 0.500:
            cmd.color("orange",   sel)
        else:
            cmd.color("red",      sel)


def draw_distance_line(atom1, atom2, name, color="red", gap=0.0):
    """
    Draw a distance line between two atoms (salt bridge / H-bond).
    atom1, atom2: PyMOL selection strings e.g. "/obj/A/LYS332/NZ"
    """
    cmd.distance(name, atom1, atom2)
    cmd.set("dash_color", color, name)
    cmd.set("dash_width", 2.5,   name)
    cmd.set("dash_radius", 0.05, name)
    cmd.set("dash_gap",    gap,  name)
    cmd.hide("labels", name)


# ─────────────────────────────────────────────────────────────
# VIEW 1: Full complex overview — cartoon + surface inset
# Mimics slide 3 of presentation
# ─────────────────────────────────────────────────────────────
def view_full_complex(pdb_file, chain_a, chain_b, anchor_resi,
                      cluster_resi_list, out_prefix):
    """
    Full complex: Chain A as surface, Chain B as cartoon.
    Anchor residue as sphere. Cluster as sticks.
    """
    setup_pymol_prism()
    cmd.load(pdb_file, "complex")
    cmd.remove("solvent")
    cmd.remove("HETATM")

    # Representations
    cmd.hide("everything", "complex")
    cmd.show("surface",  f"complex and chain {chain_a}")
    cmd.show("cartoon",  f"complex and chain {chain_b}")

    # Colors
    r, g, b = COLORS["surface_main"]
    cmd.set_color("col_surf", [r, g, b])
    cmd.color("col_surf", f"complex and chain {chain_a}")

    r, g, b = COLORS["chain_b"]
    cmd.set_color("col_cb", [r, g, b])
    cmd.color("col_cb", f"complex and chain {chain_b}")

    # Anchor residue — sphere on surface
    anchor_sel = f"complex and chain {chain_a} and resi {anchor_resi}"
    cmd.show("spheres", anchor_sel)
    cmd.color("red", anchor_sel)

    # Cluster — sticks visible through surface
    for resi in cluster_resi_list:
        sel = f"complex and chain {chain_a} and resi {resi}"
        cmd.show("sticks", sel)
        cmd.color("orange", sel)

    # Surface transparency to see cluster
    cmd.set("transparency", 0.35, f"complex and chain {chain_a}")

    # Orient and zoom
    cmd.orient("complex")
    cmd.zoom("complex", 5)

    render_figure(f"{out_prefix}_view1_full_complex.png")
    cmd.save(f"{out_prefix}_view1_full_complex.pse")
    cmd.delete("all")


# ─────────────────────────────────────────────────────────────
# VIEW 2: Binding pocket zoom — mimics slide 5 center panel
# ─────────────────────────────────────────────────────────────
def view_pocket_zoom(pdb_file, chain_a, chain_b, anchor_resi,
                     cluster_resi_list, out_prefix, zoom_radius=12):
    """
    Zoom into binding pocket. Surface + sticks for interface residues.
    Matches the zoom inset style in the presentation.
    """
    setup_pymol_prism()
    cmd.load(pdb_file, "complex")
    cmd.remove("solvent")
    cmd.remove("HETATM")

    # Show everything as surface initially
    cmd.hide("everything")
    cmd.show("surface", "complex")

    r, g, b = COLORS["surface_main"]
    cmd.set_color("col_s", [r, g, b])
    cmd.color("col_s", "complex")
    cmd.set("transparency", 0.20, "complex")

    # Cluster residues as sticks on top of surface
    for resi in cluster_resi_list:
        sel_a = f"complex and chain {chain_a} and resi {resi}"
        sel_b = f"complex and chain {chain_b} and resi {resi}"
        for sel in [sel_a, sel_b]:
            cmd.show("sticks", sel)

    # Anchor as sphere
    anchor_sel = f"complex and chain {chain_a} and resi {anchor_resi}"
    cmd.show("spheres", anchor_sel)
    cmd.color("red", anchor_sel)
    cmd.set("sphere_scale", 0.5, anchor_sel)

    # Color by element for sticks (Prism: C=tan, N=blue, O=red, S=yellow)
    cmd.util.cbay("complex")   # color by atom type, yellow carbons for ligand-like
    cmd.color("red", anchor_sel)

    # Zoom to anchor
    cmd.zoom(anchor_sel, zoom_radius)
    cmd.clip("slab", 25)       # depth clip for clean pocket view

    render_figure(f"{out_prefix}_view2_pocket_zoom.png")
    cmd.save(f"{out_prefix}_view2_pocket_zoom.pse")
    cmd.delete("all")


# ─────────────────────────────────────────────────────────────
# VIEW 3: Conservation coloring — cluster residues by cons score
# ─────────────────────────────────────────────────────────────
def view_conservation(pdb_file, chain_a, chain_b,
                      anchor_resi, cluster_resi_list,
                      conservation_dict, out_prefix):
    """
    Cluster residues colored by conservation score.
    Green = pan-coronavirus. Red = SARS-selective.
    Background protein = gray surface.
    """
    setup_pymol_prism()
    cmd.load(pdb_file, "complex")
    cmd.remove("solvent")
    cmd.remove("HETATM")

    cmd.hide("everything")
    cmd.show("surface",  "complex")
    cmd.color("gray85",  "complex")
    cmd.set("transparency", 0.50, "complex")

    # Cluster as sticks colored by conservation
    all_cluster = cluster_resi_list + [anchor_resi]
    for resi in all_cluster:
        sel_a = f"complex and chain {chain_a} and resi {resi}"
        cmd.show("sticks", sel_a)
        score = conservation_dict.get(resi, 0.5)
        if score >= 0.999:
            cmd.color("forest",  sel_a)
        elif score >= 0.800:
            cmd.color("green",   sel_a)
        elif score >= 0.500:
            cmd.color("orange",  sel_a)
        else:
            cmd.color("red",     sel_a)
        cmd.show("spheres", f"{sel_a} and name CA")
        cmd.set("sphere_scale", 0.30, f"{sel_a} and name CA")

    # Label anchor
    anchor_sel = f"complex and chain {chain_a} and resi {anchor_resi} and name CA"
    cmd.label(anchor_sel, '"%s%s" % (resn, resi)')

    cmd.zoom(f"complex and chain {chain_a} and resi {'+'.join(map(str, all_cluster))}", 10)

    render_figure(f"{out_prefix}_view3_conservation.png")
    cmd.save(f"{out_prefix}_view3_conservation.pse")
    cmd.delete("all")


# ─────────────────────────────────────────────────────────────
# VIEW 4: Salt bridge / key interaction zoom
# Mimics slide 5-6 of presentation
# ─────────────────────────────────────────────────────────────
def view_key_interaction(pdb_file, chain_a, chain_b,
                         res_a, res_b, atom_a, atom_b,
                         interaction_type, out_prefix):
    """
    Zoom into a specific salt bridge or H-bond.
    Draws dashed line between partners.
    interaction_type: 'salt_bridge' or 'hbond'
    """
    setup_pymol_prism()
    cmd.load(pdb_file, "complex")
    cmd.remove("solvent")
    cmd.remove("HETATM")

    cmd.hide("everything")
    cmd.show("cartoon", "complex")
    cmd.color("gray80", "complex")
    cmd.set("cartoon_transparency", 0.60)

    # The two interacting residues
    sel_a = f"complex and chain {chain_a} and resi {res_a}"
    sel_b = f"complex and chain {chain_b} and resi {res_b}"
    cmd.show("sticks",  f"{sel_a} or {sel_b}")
    cmd.show("spheres", f"{sel_a} or {sel_b}")
    cmd.set("sphere_scale", 0.20)
    cmd.util.cbaw(f"{sel_a} or {sel_b}")  # color by element, white C

    # Color chains distinctly
    r, g, b = COLORS["chain_a"]
    cmd.set_color("col_a", [r, g, b])
    r, g, b = COLORS["chain_b"]
    cmd.set_color("col_b2", [r, g, b])
    cmd.color("col_a",  f"{sel_a} and elem C")
    cmd.color("col_b2", f"{sel_b} and elem C")

    # Draw interaction line
    atom1 = f"/{pdb_file[:-4].split('/')[-1]}/{chain_a}/{res_a}/{atom_a}"
    atom2 = f"/{pdb_file[:-4].split('/')[-1]}/{chain_b}/{res_b}/{atom_b}"
    color = "red" if interaction_type == "salt_bridge" else "blue"
    draw_distance_line(atom1, atom2, "interaction", color=color, gap=0.2)
    cmd.show("dashes", "interaction")

    # Labels
    cmd.label(f"{sel_a} and name CA", '"%s%s" % (resn, resi)')
    cmd.label(f"{sel_b} and name CA", '"%s%s" % (resn, resi)')

    # Zoom
    center_sel = f"{sel_a} or {sel_b}"
    cmd.zoom(center_sel, 8)
    cmd.clip("slab", 15)

    render_figure(f"{out_prefix}_view4_{interaction_type}.png")
    cmd.save(f"{out_prefix}_view4_{interaction_type}.pse")
    cmd.delete("all")


# ─────────────────────────────────────────────────────────────
# VIEW 5: Pharmacophore cluster spheres
# Shows the abstract pharmacophore hypothesis in 3D
# ─────────────────────────────────────────────────────────────
def view_pharmacophore(pdb_file, chain_a, chain_b,
                       pharmacophore_dict, out_prefix):
    """
    Show pharmacophore features as colored spheres at residue centroids.
    pharmacophore_dict: {resi: feature_type}
    feature_type: 'anchor' | 'aromatic' | 'hbond_donor' |
                  'hbond_acceptor' | 'hydrophobic' | 'charged_pos' |
                  'charged_neg'
    """
    feature_colors = {
        "anchor":        "red",
        "aromatic":      "orange",
        "hbond_donor":   "blue",
        "hbond_acceptor": "cyan",
        "hydrophobic":   "yellow",
        "charged_pos":   "marine",
        "charged_neg":   "salmon",
    }

    setup_pymol_prism()
    cmd.load(pdb_file, "complex")
    cmd.remove("solvent")
    cmd.remove("HETATM")

    cmd.hide("everything")
    cmd.show("surface", "complex")
    cmd.color("gray90", "complex")
    cmd.set("transparency", 0.55, "complex")
    cmd.show("cartoon",  "complex")
    cmd.color("gray70",  "complex")
    cmd.set("cartoon_transparency", 0.70)

    # Draw pharmacophore spheres at Cα positions
    for resi, feature in pharmacophore_dict.items():
        sel = f"complex and chain {chain_a} and resi {resi} and name CA"
        cmd.show("spheres", sel)
        cmd.set("sphere_scale", 0.65, sel)
        color = feature_colors.get(feature, "white")
        cmd.color(color, sel)
        cmd.label(sel, '"%s" % resn')

    resi_list = list(pharmacophore_dict.keys())
    zoom_sel = f"complex and chain {chain_a} and resi {'+'.join(map(str, resi_list))}"
    cmd.zoom(zoom_sel, 12)

    render_figure(f"{out_prefix}_view5_pharmacophore.png")
    cmd.save(f"{out_prefix}_view5_pharmacophore.pse")
    cmd.delete("all")


# ─────────────────────────────────────────────────────────────
# VIEW 6: Multi-structure overlay (crystal + AF3)
# Mimics the structural overlay views in the pipeline
# ─────────────────────────────────────────────────────────────
def view_structural_overlay(pdb_files_dict, chain_a, anchor_resi,
                             out_prefix):
    """
    Overlay multiple structures of same interface.
    pdb_files_dict: {'7BV2': path, '6NUR': path, 'AF3': path}
    Colors each structure distinctly.
    """
    overlay_colors = ["blue", "forest", "red", "orange", "purple"]

    setup_pymol_prism()

    for i, (name, path) in enumerate(pdb_files_dict.items()):
        cmd.load(path, name)
        cmd.hide("everything", name)
        cmd.show("cartoon",    name)
        color = overlay_colors[i % len(overlay_colors)]
        cmd.color(color, name)
        cmd.set("cartoon_transparency", 0.30 if i > 0 else 0.0, name)

    # Align all to first
    names = list(pdb_files_dict.keys())
    for name in names[1:]:
        cmd.align(name, names[0])

    # Highlight anchor in all structures
    for name in names:
        anchor_sel = f"{name} and chain {chain_a} and resi {anchor_resi}"
        cmd.show("spheres", anchor_sel)
        cmd.set("sphere_scale", 0.4, anchor_sel)

    cmd.orient(names[0])
    cmd.zoom(f"{names[0]} and chain {chain_a} and resi {anchor_resi}", 15)

    render_figure(f"{out_prefix}_view6_overlay.png")
    cmd.save(f"{out_prefix}_view6_overlay.pse")
    cmd.delete("all")


# ─────────────────────────────────────────────────────────────
# VIEW 7: Control comparison panel
# Positive control (known hotspot) vs Negative control (surface residue)
# This is what was MISSING from the JB12/JB16 presentation
# ─────────────────────────────────────────────────────────────
def view_controls(pdb_file, chain_a, chain_b,
                  positive_control_resi,   # Known critical residue
                  negative_control_resi,   # Known non-critical residue
                  out_prefix):
    """
    Side-by-side control comparison.
    Left: positive control residue (known hotspot — should be buried, large BSA)
    Right: negative control residue (surface residue — should be exposed, low BSA)

    This proves the pipeline correctly distinguishes interface from surface.
    """
    for ctrl_type, resi, out_label in [
        ("positive", positive_control_resi, "pos_ctrl"),
        ("negative", negative_control_resi, "neg_ctrl"),
    ]:
        setup_pymol_prism()
        cmd.load(pdb_file, "complex")
        cmd.remove("solvent")
        cmd.remove("HETATM")

        cmd.hide("everything")
        cmd.show("surface",  f"complex and chain {chain_a}")
        cmd.color("gray80",  f"complex and chain {chain_a}")
        cmd.set("transparency", 0.30)

        cmd.show("cartoon",  f"complex and chain {chain_b}")
        r, g, b = COLORS["chain_b"]
        cmd.set_color("ctrl_cb", [r, g, b])
        cmd.color("ctrl_cb",  f"complex and chain {chain_b}")

        # Highlight control residue
        ctrl_sel = f"complex and chain {chain_a} and resi {resi}"
        cmd.show("sticks",  ctrl_sel)
        cmd.show("spheres", f"{ctrl_sel} and name CA")
        cmd.set("sphere_scale", 0.50, f"{ctrl_sel} and name CA")

        if ctrl_type == "positive":
            cmd.color("red", ctrl_sel)   # positive = red, should be buried
        else:
            cmd.color("blue", ctrl_sel)  # negative = blue, should be exposed

        cmd.label(f"{ctrl_sel} and name CA", '"%s%s" % (resn, resi)')
        cmd.zoom(ctrl_sel, 12)

        render_figure(f"{out_prefix}_view7_{out_label}.png")
        cmd.delete("all")


# ─────────────────────────────────────────────────────────────
# GENERATE ALL VIEWS — main function
# ─────────────────────────────────────────────────────────────
def generate_all_views(
    complex_name,
    anchor_resi,
    pdb_file,
    chain_a,
    chain_b,
    cluster_resi_list,
    conservation_dict,
    pharmacophore_dict,
    salt_bridge_pairs,        # list of (res_a, atom_a, res_b, atom_b)
    positive_control_resi,
    negative_control_resi,
    pdb_files_overlay,        # dict for multi-structure overlay
    output_dir
):
    """
    Generate all 7 standard views for one interface.
    All figures saved to output_dir.
    """
    os.makedirs(output_dir, exist_ok=True)
    prefix = os.path.join(output_dir, complex_name)

    print(f"\n{'='*60}")
    print(f"Generating 3D figures for: {complex_name}")
    print(f"Anchor: {anchor_resi} | Cluster: {cluster_resi_list}")
    print(f"Output: {output_dir}")
    print(f"{'='*60}\n")

    print("View 1: Full complex overview...")
    view_full_complex(pdb_file, chain_a, chain_b,
                      anchor_resi, cluster_resi_list, prefix)

    print("View 2: Pocket zoom...")
    view_pocket_zoom(pdb_file, chain_a, chain_b,
                     anchor_resi, cluster_resi_list, prefix)

    print("View 3: Conservation coloring...")
    view_conservation(pdb_file, chain_a, chain_b,
                      anchor_resi, cluster_resi_list,
                      conservation_dict, prefix)

    print("View 4: Key interaction (salt bridge / H-bond)...")
    for i, (ra, aa, rb, ab) in enumerate(salt_bridge_pairs):
        view_key_interaction(pdb_file, chain_a, chain_b,
                             ra, rb, aa, ab,
                             "salt_bridge", f"{prefix}_sb{i+1}")

    print("View 5: Pharmacophore hypothesis...")
    view_pharmacophore(pdb_file, chain_a, chain_b,
                       pharmacophore_dict, prefix)

    print("View 6: Multi-structure overlay...")
    if pdb_files_overlay:
        view_structural_overlay(pdb_files_overlay, chain_a,
                                anchor_resi, prefix)

    print("View 7: Control comparison...")
    view_controls(pdb_file, chain_a, chain_b,
                  positive_control_resi,
                  negative_control_resi,
                  prefix)

    print(f"\nAll views complete. Files in: {output_dir}")


# ─────────────────────────────────────────────────────────────
# EXAMPLE USAGE — NSP12-NSP7
# Run: pymol -c pymol_cgcp_visuals.py
# ─────────────────────────────────────────────────────────────
if __name__ == "pymol":

    BASE = os.path.expanduser("~/projects/rtc-pan-coronavirus")

    generate_all_views(
        complex_name       = "NSP12-NSP7",
        anchor_resi        = 440,           # PHE440 — primary anchor
        pdb_file           = f"{BASE}/00-reference/known_interfaces/NSP12-NSP7/7BV2_NSP12-NSP7.pdb",
        chain_a            = "A",           # NSP12
        chain_b            = "C",           # NSP7
        cluster_resi_list  = [412, 413, 415, 420, 442, 843],
        conservation_dict  = {
            440: 1.000,   # PHE440 — pan-coronavirus
            412: 1.000,   # PRO412
            413: 1.000,   # GLY413
            415: 1.000,   # PHE415
            420: 1.000,   # TYR420
            442: 1.000,   # PHE442
            843: 1.000,   # PHE843
            431: 1.000,   # GLU431 — salt bridge
        },
        pharmacophore_dict = {
            440: "anchor",         # PHE440 — aromatic anchor
            412: "hydrophobic",    # PRO412
            415: "aromatic",       # PHE415
            420: "hbond_donor",    # TYR420
            442: "aromatic",       # PHE442
            431: "charged_neg",    # GLU431 — salt bridge
        },
        salt_bridge_pairs  = [
            (431, "OE1", 1, "NZ"),  # GLU431(NSP12) — LYS1(NSP7)
        ],
        positive_control_resi = 440,   # PHE440 — known anchor, should be buried
        negative_control_resi = 100,   # Surface residue — should show low contacts
        pdb_files_overlay  = {
            "7BV2": f"{BASE}/00-reference/known_interfaces/NSP12-NSP7/7BV2_NSP12-NSP7.pdb",
            "6NUR": f"{BASE}/00-reference/known_interfaces/NSP12-NSP7/6NUR_NSP12-NSP7.pdb",
            "7C2K": f"{BASE}/00-reference/known_interfaces/NSP12-NSP7/7C2K_NSP12-NSP7.pdb",
        },
        output_dir = f"{BASE}/CGCP/02-deep-dive/NSP12-NSP7/step-01-structure/visuals"
    )

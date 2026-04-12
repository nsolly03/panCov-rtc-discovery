#!/usr/bin/env python3
"""
CGCP Docking Box Verification Script v2
=========================================
Strategy:
  1. Read step-06 assessment TSV → get ANCHOR + INCLUDE residues + their
     sequence positions (the `position` column = PDB residue seq number)
  2. Look up Cα coordinates from interface-specific PDB files
  3. Compute true bounding box of all pharmacophore residues + 10Å buffer
  4. Compare to current v2 docking box
  5. Report PASS / FAIL + corrected coordinates ready for VFUparr all.ctrl

Usage:
    conda activate rtc-discovery
    cd ~/projects/rtc-pan-coronavirus
    python3 CGCP/scripts/phase-2/verify_docking_boxes.py

Outputs:
    CGCP/02-deep-dive/box-verification/verify_docking_boxes_report.tsv
    CGCP/02-deep-dive/box-verification/corrected_boxes_vfuparr.txt
    CGCP/02-deep-dive/box-verification/Fig_BoxVerification.png

Author : CGCP Pipeline | Olivier Nsekuye | ULiège GIGA-VIN
Date   : 2026-04-12
"""

import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────
PROJECT  = os.path.expanduser("~/projects/rtc-pan-coronavirus")
DEEP     = os.path.join(PROJECT, "CGCP", "02-deep-dive")
PDB_DIR  = os.path.join(PROJECT, "00-reference", "pdb_structures")
OUT_DIR  = os.path.join(DEEP, "box-verification")
os.makedirs(OUT_DIR, exist_ok=True)

BUFFER   = 10.0   # Å added to each side of bounding box
DECISIONS_KEEP = {"ANCHOR", "ANCHOR_PRIMARY", "ANCHOR_SECONDARY", "INCLUDE"}

# ─────────────────────────────────────────────────────────────────────────────
# INTERFACE REGISTRY
# ─────────────────────────────────────────────────────────────────────────────
INTERFACES = {
    "NSP12-NSP7": {
        "tsv"       : "NSP12-NSP7/step-06-assessment/integrated_assessment_NSP12-NSP7.tsv",
        "pdb"       : "03-virtual-screening/NSP12-NSP7_3/receptor_NSP12-NSP7_3.pdb",
        "chain_map" : {"NSP12": "A", "NSP7": "C"},
        "v2_box"    : {"center": (100.015,  89.131, 106.015),
                       "size"  : ( 32.235,  58.390,  54.642)},
    },
    "NSP12-NSP8": {
        "tsv"       : "NSP12-NSP8/step-06-assessment/integrated_assessment_NSP12-NSP8.tsv",
        "pdb"       : "03-virtual-screening/NSP12-NSP8_4/receptor_NSP12-NSP8_4.pdb",
        "chain_map" : {"NSP12": "A", "NSP8": "B"},
        "v2_box"    : {"center": ( 98.620, 119.740, 113.310),
                       "size"  : ( 65.000,  49.000,  44.000)},
    },
    "NSP9-NSP12": {
        "tsv"       : "NSP9-NSP12/step-06-assessment/integrated_assessment_NSP9-NSP12.tsv",
        "pdb"       : "03-virtual-screening/NSP9-NSP12_5/receptor_NSP9-NSP12_5.pdb",
        "chain_map" : {"NSP12": "A", "NSP9": "G"},
        "v2_box"    : {"center": (133.250, 167.400, 176.140),
                       "size"  : ( 46.000,  38.000,  46.000)},
    },
    "NSP10-NSP16": {
        "tsv"       : "NSP10-NSP16/step-06-assessment/integrated_assessment_NSP10-NSP16.tsv",
        "pdb"       : "03-virtual-screening/NSP10-NSP16_2/receptor_NSP10-NSP16_2.pdb",
        # TSV position column uses polyprotein absolute numbers — matches PDB directly
        "chain_map" : {"NSP16": "A", "NSP10": "B"},
        "v2_box"    : {"center": ( 73.510,  22.270,  14.860),
                       "size"  : ( 32.000,  44.000,  39.000)},
    },
    "NSP7-NSP8": {
        "tsv"       : "NSP7-NSP8/step-06-assessment/integrated_assessment_NSP7-NSP8.tsv",
        "pdb"       : "03-virtual-screening/NSP7-NSP8_6/receptor_NSP7-NSP8_ModeB_AF3_6.pdb",
        "chain_map" : {"NSP8": "A", "NSP7": "B"},
        "v2_box"    : {"center": ( -7.980, -11.619, -11.870),
                       "size"  : ( 47.300,  39.900,  37.000)},
    },
    "NSP10-NSP14": {
        "tsv"       : "NSP10-NSP14/step-06-assessment/integrated_assessment_NSP10-NSP14.tsv",
        "pdb"       : "03-virtual-screening/NSP10-NSP14_2/receptor_NSP10-NSP14_2.pdb",
        "chain_map" : {"NSP10": "B", "NSP14": "A"},
        "v2_box"    : {"center": ( -3.700,   7.460, -25.650),
                       "size"  : ( 39.000,  42.000,  57.000)},
    },
    "NSP13-Helicase": {
        "tsv"       : "NSP13-Helicase/step-06-assessment/integrated_assessment_NSP13-Helicase.tsv",
        "pdb"       : "03-virtual-screening/NSP13-Helicase_7/receptor_NSP13-Helicase_7.pdb",
        "chain_map" : {"NSP13a": "A", "NSP13b": "E"},
        "v2_box"    : {"center": (-30.151,  14.648,  -9.240),
                       "size"  : ( 46.000,  46.000,  36.000)},
    },
    "NSP12-NSP13": {
        "tsv"       : "NSP12-NSP13/step-06-assessment/integrated_assessment_NSP12-NSP13.tsv",
        "pdb"       : "03-virtual-screening/NSP12-NSP13_8/receptor_NSP12-NSP13_8.pdb",
        "chain_map" : {"NSP12": "A", "NSP13": "E"},
        "v2_box"    : {"center": (198.500, 194.000, 159.260),
                       "size"  : ( 24.000,  24.000,  24.000)},
    },
    "NSP15": {
        "tsv"       : "NSP15/step-06-assessment/integrated_assessment_NSP15.tsv",
        "pdb"       : "03-virtual-screening/NSP15_9/receptor_NSP15_9.pdb",
        "chain_map" : {"NSP15a": "A", "NSP15b": "B"},
        "v2_box"    : {"center": (-70.407,  44.044, -15.772),
                       "size"  : ( 52.000,  52.000,  50.000)},
    },
}

# ─────────────────────────────────────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def load_tsv_residues(tsv_path):
    if not os.path.exists(tsv_path):
        return None, f"TSV not found: {tsv_path}"
    residues = []
    with open(tsv_path) as f:
        header = None
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if header is None:
                header = [c.strip() for c in cols]
                continue
            row = dict(zip(header, cols))
            dec = row.get("decision", "").strip().upper()
            if dec in DECISIONS_KEEP:
                residues.append({
                    "residue" : row.get("residue", "").strip(),
                    "chain"   : row.get("chain", "").strip(),
                    "position": row.get("position", "").strip(),
                    "aa"      : row.get("aa", "").strip(),
                    "feature" : row.get("primary_feature", "").strip(),
                    "decision": dec,
                })
    if not residues:
        return None, "No ANCHOR/INCLUDE residues found"
    return residues, None


def load_pdb_ca(pdb_path):
    coords = {}
    if not os.path.exists(pdb_path):
        return coords, f"PDB not found: {pdb_path}"
    with open(pdb_path) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            if line[12:16].strip() != "CA":
                continue
            chain = line[21].strip()
            try:
                res_seq = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            coords[(chain, res_seq)] = (x, y, z)
    return coords, None


def match_residues(residues, pdb_coords, chain_map):
    found, missing = {}, []
    for r in residues:
        try:
            pos = int(r["position"])
        except ValueError:
            missing.append(r["residue"]); continue

        coord = None
        pdb_chain = chain_map.get(r["chain"])
        if pdb_chain:
            coord = pdb_coords.get((pdb_chain, pos))
        if coord is None:
            for pc in chain_map.values():
                coord = pdb_coords.get((pc, pos))
                if coord: break
        if coord is None:
            for (ch, seq), xyz in pdb_coords.items():
                if seq == pos:
                    coord = xyz; break
        if coord:
            found[r["residue"]] = coord
        else:
            missing.append(r["residue"])
    return found, missing


def bounding_box(coords_dict, buffer=BUFFER):
    pts    = np.array(list(coords_dict.values()))
    mn, mx = pts.min(axis=0), pts.max(axis=0)
    center = (mn + mx) / 2.0
    size   = (mx - mn) + 2 * buffer
    return tuple(center.tolist()), tuple(size.tolist())


def check_coverage(v2_box, req_center, req_size):
    vc, vs = np.array(v2_box["center"]), np.array(v2_box["size"])
    rc, rs = np.array(req_center),       np.array(req_size)
    v_min, v_max = vc - vs/2, vc + vs/2
    r_min, r_max = rc - rs/2, rc + rs/2
    issues = []
    for i, ax in enumerate(["X","Y","Z"]):
        if v_min[i] > r_min[i]:
            issues.append(f"{ax}-min: v2={v_min[i]:.2f} > req={r_min[i]:.2f} (under {v_min[i]-r_min[i]:.2f}Å)")
        if v_max[i] < r_max[i]:
            issues.append(f"{ax}-max: v2={v_max[i]:.2f} < req={r_max[i]:.2f} (short {r_max[i]-v_max[i]:.2f}Å)")
    return issues


def round_box(center, size):
    c = tuple(round(v, 3) for v in center)
    s = tuple(float(np.ceil(v / 2) * 2) for v in size)
    return c, s


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    print("\n" + "="*72)
    print("  CGCP DOCKING BOX VERIFICATION  v2  —  All 9 Interfaces")
    print("="*72)

    results = []

    for iface, cfg in INTERFACES.items():
        print(f"\n{'─'*60}\n  {iface}\n{'─'*60}")

        # 1. TSV
        tsv_path = os.path.join(DEEP, cfg["tsv"])
        residues, err = load_tsv_residues(tsv_path)
        if err:
            print(f"  ⚠  {err}")
            results.append({"interface":iface,"status":"TSV_ERROR","issues":[err]}); continue
        print(f"  Pharmacophore residues (ANCHOR+INCLUDE): {len(residues)}")

        # 2. PDB
        pdb_path = os.path.join(PROJECT, cfg["pdb"])
        pdb_coords, err = load_pdb_ca(pdb_path)
        if err:
            print(f"  ⚠  {err}")
            results.append({"interface":iface,"status":"PDB_ERROR","issues":[err]}); continue
        print(f"  PDB CA atoms loaded: {len(pdb_coords)}  ({os.path.basename(pdb_path)})")

        # 3. Match
        found, missing = match_residues(residues, pdb_coords, cfg["chain_map"])
        print(f"  Residues matched: {len(found)}/{len(residues)}")
        if missing:
            print(f"  Missing: {', '.join(missing[:8])}" + (" ..." if len(missing)>8 else ""))

        if len(found) < 2:
            print(f"  ⚠  Too few coordinates (<2) to compute box")
            results.append({"interface":iface,"status":"COORD_ERROR",
                             "issues":["<2 residues matched"]}); continue

        # 4. Required box
        req_center, req_size = bounding_box(found)
        print(f"  Required box:  center=({req_center[0]:.2f}, {req_center[1]:.2f}, {req_center[2]:.2f})"
              f"  size={req_size[0]:.1f}×{req_size[1]:.1f}×{req_size[2]:.1f}Å")

        # 5. v2 box
        v2 = cfg["v2_box"]
        print(f"  v2 box:        center=({v2['center'][0]:.2f}, {v2['center'][1]:.2f}, {v2['center'][2]:.2f})"
              f"  size={v2['size'][0]:.1f}×{v2['size'][1]:.1f}×{v2['size'][2]:.1f}Å")

        # 6. Check
        issues      = check_coverage(v2, req_center, req_size)
        center_disp = float(np.linalg.norm(np.array(v2["center"]) - np.array(req_center)))
        sug_c, sug_s = round_box(req_center, req_size)
        v2_vol  = v2["size"][0]  * v2["size"][1]  * v2["size"][2]
        req_vol = req_size[0]    * req_size[1]    * req_size[2]

        print(f"  Centre displacement: {center_disp:.2f}Å  |  "
              f"Volume v2/req: {v2_vol/1000:.1f}k / {req_vol/1000:.1f}k Å³")

        if not issues:
            status = "PASS"
            print(f"  ✅ PASS — v2 box fully covers all pharmacophore elements")
        else:
            status = "FAIL"
            print(f"  ❌ FAIL — {len(issues)} issue(s):")
            for iss in issues:
                print(f"     • {iss}")
            print(f"  Suggested: center=({sug_c[0]:.3f},{sug_c[1]:.3f},{sug_c[2]:.3f})"
                  f"  size={sug_s[0]:.1f}×{sug_s[1]:.1f}×{sug_s[2]:.1f}Å")

        results.append({"interface":iface,"status":status,"issues":issues,
                        "n_residues":len(residues),"n_matched":len(found),
                        "n_missing":len(missing),
                        "req_center":req_center,"req_size":req_size,
                        "v2_center":v2["center"],"v2_size":v2["size"],
                        "sug_center":sug_c,"sug_size":sug_s,
                        "center_disp":center_disp,
                        "v2_vol":v2_vol,"req_vol":req_vol})

    # ── Summary ───────────────────────────────────────────────────────────
    print("\n" + "="*72 + "\n  SUMMARY\n" + "="*72)
    print(f"  {'Interface':<20} {'Status':<8} {'Matched':>8} {'Δcenter':>9}  Issues")
    print(f"  {'─'*20} {'─'*8} {'─'*8} {'─'*9}  {'─'*30}")
    for r in results:
        sym  = "✅" if r["status"]=="PASS" else ("❌" if r["status"]=="FAIL" else "⚠")
        mstr = f"{r.get('n_matched',0)}/{r.get('n_residues',0)}"
        dstr = f"{r.get('center_disp',0):.1f}Å"
        print(f"  {r['interface']:<20} {sym} {r['status']:<6} {mstr:>8} {dstr:>9}  "
              f"{len(r.get('issues',[]))} issue(s)")

    pass_n = sum(1 for r in results if r["status"]=="PASS")
    fail_n = sum(1 for r in results if r["status"]=="FAIL")
    print(f"\n  Result: {pass_n} PASS  |  {fail_n} FAIL  |  {len(results)-pass_n-fail_n} ERROR")

    # ── VFUparr coordinate table ──────────────────────────────────────────
    print("\n" + "="*72)
    print("  FINAL COORDINATES FOR VFUparr all.ctrl")
    print("="*72)
    print(f"  {'Interface':<20} {'CENTER-X':>10} {'CENTER-Y':>10} {'CENTER-Z':>10}"
          f" {'SIZE-X':>8} {'SIZE-Y':>8} {'SIZE-Z':>8}  Action")
    print(f"  {'─'*20} {'─'*10} {'─'*10} {'─'*10} {'─'*8} {'─'*8} {'─'*8}  {'─'*12}")

    vfuparr_blocks = []
    for r in results:
        if r["status"] not in ("PASS","FAIL"): continue
        c = r["v2_center"] if r["status"]=="PASS" else r["sug_center"]
        s = r["v2_size"]   if r["status"]=="PASS" else r["sug_size"]
        action = "USE v2" if r["status"]=="PASS" else "⚠ CORRECTED"
        print(f"  {r['interface']:<20} {c[0]:>10.3f} {c[1]:>10.3f} {c[2]:>10.3f}"
              f" {s[0]:>8.1f} {s[1]:>8.1f} {s[2]:>8.1f}  {action}")
        vfuparr_blocks.append(
            f"# {r['interface']}  [{action}]\n"
            f"CENTER-X={c[0]:.3f}\nCENTER-Y={c[1]:.3f}\nCENTER-Z={c[2]:.3f}\n"
            f"SIZE-X={s[0]:.1f}\nSIZE-Y={s[1]:.1f}\nSIZE-Z={s[2]:.1f}\n"
        )

    # ── Write files ───────────────────────────────────────────────────────
    tsv_out = os.path.join(OUT_DIR, "verify_docking_boxes_report.tsv")
    with open(tsv_out, "w") as f:
        f.write("Interface\tStatus\tN_residues\tN_matched\tN_missing\tCenter_disp_A\t"
                "Req_CX\tReq_CY\tReq_CZ\tReq_SX\tReq_SY\tReq_SZ\t"
                "V2_CX\tV2_CY\tV2_CZ\tV2_SX\tV2_SY\tV2_SZ\t"
                "Sug_CX\tSug_CY\tSug_CZ\tSug_SX\tSug_SY\tSug_SZ\tIssues\n")
        for r in results:
            rc=r.get("req_center",(0,0,0)); rs=r.get("req_size",(0,0,0))
            vc=r.get("v2_center",(0,0,0));  vs=r.get("v2_size",(0,0,0))
            sc=r.get("sug_center",(0,0,0)); ss=r.get("sug_size",(0,0,0))
            f.write(f"{r['interface']}\t{r['status']}\t"
                    f"{r.get('n_residues','')}\t{r.get('n_matched','')}\t{r.get('n_missing','')}\t"
                    f"{r.get('center_disp',0):.3f}\t"
                    f"{rc[0]:.3f}\t{rc[1]:.3f}\t{rc[2]:.3f}\t"
                    f"{rs[0]:.1f}\t{rs[1]:.1f}\t{rs[2]:.1f}\t"
                    f"{vc[0]:.3f}\t{vc[1]:.3f}\t{vc[2]:.3f}\t"
                    f"{vs[0]:.1f}\t{vs[1]:.1f}\t{vs[2]:.1f}\t"
                    f"{sc[0]:.3f}\t{sc[1]:.3f}\t{sc[2]:.3f}\t"
                    f"{ss[0]:.1f}\t{ss[1]:.1f}\t{ss[2]:.1f}\t"
                    f"{'|'.join(r.get('issues',[]))}\n")

    vf_out = os.path.join(OUT_DIR, "corrected_boxes_vfuparr.txt")
    with open(vf_out, "w") as f:
        f.write("# CGCP VFUparr Docking Box Coordinates\n"
                "# Generated by verify_docking_boxes.py v2  |  2026-04-12\n\n")
        f.write("\n".join(vfuparr_blocks))

    print(f"\n  TSV    → {tsv_out}")
    print(f"  VFUparr→ {vf_out}")

    # ── Figure ────────────────────────────────────────────────────────────
    fig_data = [r for r in results if r.get("v2_vol") and r.get("req_vol")]
    if fig_data:
        labels = [r["interface"].replace("Helicase","Hel") for r in fig_data]
        x, w   = np.arange(len(labels)), 0.35
        c_v2   = "#4472C4"
        col_req= ["#FF0000" if r["status"]=="FAIL" else "#70AD47" for r in fig_data]

        fig, (ax1,ax2) = plt.subplots(1,2,figsize=(14,5),facecolor="white")
        fig.patch.set_facecolor("white")

        ax1.bar(x-w/2, [r["v2_vol"]/1000  for r in fig_data], w, color=c_v2,  zorder=3)
        ax1.bar(x+w/2, [r["req_vol"]/1000 for r in fig_data], w, color=col_req,zorder=3)
        ax1.set_xticks(x); ax1.set_xticklabels(labels,rotation=45,ha="right",fontsize=8)
        ax1.set_ylabel("Box volume (×10³ Å³)",fontsize=10)
        ax1.set_title("A   Box Volume: v2 vs Required",fontsize=11,fontweight="bold",loc="left")
        ax1.spines[["top","right"]].set_visible(False); ax1.tick_params(direction="out")
        ax1.legend(handles=[mpatches.Patch(color=c_v2,label="v2 box"),
                             mpatches.Patch(color="#70AD47",label="Required — PASS"),
                             mpatches.Patch(color="#FF0000",label="Required — FAIL")],
                   fontsize=8, frameon=False)

        disps  = [r.get("center_disp",0) for r in fig_data]
        colors = ["#FF0000" if r["status"]=="FAIL" else "#70AD47" for r in fig_data]
        bars   = ax2.bar(x, disps, color=colors, zorder=3, width=0.55)
        ax2.axhline(5, color="#808080", ls="--", lw=0.9, label="5Å threshold")
        ax2.set_xticks(x); ax2.set_xticklabels(labels,rotation=45,ha="right",fontsize=8)
        ax2.set_ylabel("Centre displacement (Å)",fontsize=10)
        ax2.set_title("B   Box Centre: v2 vs Pharmacophore Centroid",fontsize=11,
                      fontweight="bold",loc="left")
        ax2.spines[["top","right"]].set_visible(False); ax2.tick_params(direction="out")
        ax2.legend(fontsize=8, frameon=False)
        for bar,val in zip(bars,disps):
            ax2.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.3,
                     f"{val:.1f}", ha="center", va="bottom", fontsize=7)

        plt.tight_layout()
        fig_path = os.path.join(OUT_DIR, "Fig_BoxVerification.png")
        plt.savefig(fig_path, dpi=150, bbox_inches="tight", facecolor="white")
        plt.close()
        print(f"  Figure → {fig_path}")

    sys.exit(1 if fail_n > 0 else 0)


if __name__ == "__main__":
    main()

# CGCP Step 1 Structural Audit — NSP10-NSP16

**Result: 7/7 PASS**

**Note:** Polyprotein numbering — NSP16=ChainA(6798-7096), NSP10=ChainB(4271-4386)

| Check | Name | Result | Detail |
|-------|------|--------|--------|
| C1 | Chains A (NSP16) + B (NSP10) present | ✅ | Chains: ['A', 'B'] |
| C2 | Chain sizes reasonable | ✅ | NSP16=299 res, NSP10=116 res |
| C3 | Anchor adjacent to NSP16 | ✅ | LYS93(B4346) min dist to NSP16: 5.15Å |
| C4 | Backbone complete (Chain A/NSP16) | ✅ | 0 missing backbone atoms |
| C5 | LYS93 (B4346) present and correct | ✅ | B4346=LYS at (74.536,12.890,9.624) |
| C6 | Positive ctrl: LYS93 near NSP16 | ✅ | 24 NSP16 residues within 15Å |
| C7 | Neg ctrl: LYS93 not at Zn site | ✅ | LYS93 to Zn: 22.53Å |

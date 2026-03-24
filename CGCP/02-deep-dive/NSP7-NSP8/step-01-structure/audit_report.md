# CGCP Step 1 Structural Audit — NSP7-NSP8

**Result: 7/7 PASS**

**Mode B (AlphaFold3) used — Mode A not druggable**

| Check | Name | Result | Detail |
|-------|------|--------|--------|
| C1 | Chains A (NSP8) + B (NSP7) present | ✅ | Chains: ['A', 'B'] |
| C2 | Chain sizes reasonable | ✅ | NSP8=198 res, NSP7=83 res |
| C3 | PHE92 adjacent to NSP7 | ✅ | PHE92 min dist to NSP7: 6.78Å |
| C4 | Backbone complete (Chain A/NSP8) | ✅ | 0 missing backbone atoms |
| C5 | PHE92 present and correct | ✅ | Res 92=PHE at (-3.107,-13.031,-12.646) |
| C6 | Positive ctrl: PHE92 near NSP7 | ✅ | 36 NSP7 residues within 15Å |
| C7 | Neg ctrl: PHE92 far from NSP8 C-term | ✅ | PHE92 to NSP8-191: 28.65Å |

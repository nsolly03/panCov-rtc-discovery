# CGCP Step 1 Structural Audit — NSP12-NSP8

**Result: 7/7 PASS**

| Check | Name | Result | Detail |
|-------|------|--------|--------|
| C1 | Chains A+B present | ✅ PASS | Chains found: ['A', 'B'] |
| C2 | Chain sizes reasonable | ✅ PASS | NSP12=834 res, NSP8=114 res |
| C3 | Anchor adjacent to NSP8 | ✅ PASS | LYS332 min dist to NSP8: 7.83A |
| C4 | Backbone complete (Chain A) | ✅ PASS | 0 missing backbone atoms |
| C5 | LYS332 present and correct | ✅ PASS | Res 332 = LYS at (95.339, 131.523, 109.793) |
| C6 | Positive ctrl: LYS332 near NSP8 | ✅ PASS | 29 NSP8 residues within 15A of LYS332 |
| C7 | Neg ctrl: LYS332 far from catalytic triad | ✅ PASS | res618=48.3A | res759=42.8A | res760=43.7A |

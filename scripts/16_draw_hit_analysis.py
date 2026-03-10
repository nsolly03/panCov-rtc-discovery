from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image, ImageDraw, ImageFont
import io
from pathlib import Path

Path('figures/script16').mkdir(parents=True, exist_ok=True)

data = [
    ('351017',   'O=C1c2ccccc2[C@]2(O)Nc3[nH]c(=O)[nH]c(=O)c3[C@]12O',    (-8.153,-8.132,-7.759), '#1 Lead'),
    ('13633807', 'O=C1NC(=O)C2=Nn3c(nc4ccccc4c3=O)C[C@]2(O)N1',            (-8.081,-8.173,-7.218), '#2'),
    ('13633805', 'O=C1NC(=O)C2=Nn3c(nc4ccccc4c3=O)C[C@@]2(O)N1',           (-8.087,-8.198,-7.168), '#3'),
    ('351016',   'O=C1c2ccccc2[C@@]2(O)Nc3[nH]c(=O)[nH]c(=O)c3[C@]12O',   (-7.471,-8.407,-7.439), '#4'),
    ('5024943',  'CC1(C)C[C@@]2(CCO1)NNC(=O)[C@@H]1C(=O)NNC(=O)[C@@H]12', (-8.225,-7.261,-7.790), '#5'),
    ('351018',   'O=C1c2ccccc2[C@@]2(O)Nc3[nH]c(=O)[nH]c(=O)c3[C@@]12O',  (-7.517,-8.040,-7.553), '#6'),
    ('3169307',  'NC(NC(=O)c1ccccc1)=C1C(=O)NC(=O)NC1=O',                  (-7.346,-7.858,-7.115), '#7'),
    ('5024944',  'CC1(C)C[C@]2(CCO1)NNC(=O)[C@@H]1C(=O)NNC(=O)[C@H]12',   (-7.633,-7.479,-7.185), '#8'),
    ('351019',   'O=C1c2ccccc2[C@]2(O)Nc3[nH]c(=O)[nH]c(=O)c3[C@@]12O',   (-7.168,-7.491,-7.542), '#9'),
    ('13662104', 'CC1(C)C[C@]2(CCO1)NNC(=O)[C@H]1C(=O)NNC(=O)[C@@H]12',   (-7.545,-7.492,-7.091), '#10'),
    ('5024945',  'CC1(C)C[C@@]2(CCO1)NNC(=O)[C@@H]1C(=O)NNC(=O)[C@H]12',  (-7.185,-7.275,-7.262), '#11'),
]

COLS     = 4
MOL_W    = 400
MOL_H    = 290
LABEL_H  = 105
PAD      = 18
HEADER_H = 110
FOOTER_H = 55
ROWS     = -(-len(data) // COLS)

# Light palette
BG        = (248, 249, 252)   # near-white page
CARD_BG   = (255, 255, 255)   # pure white cards
CARD_BG_A = (255, 252, 240)   # warm tint — Scaffold A
CARD_BG_B = (240, 252, 245)   # green tint — Scaffold B
CARD_BG_C = (240, 246, 255)   # blue tint — Scaffold C
LEAD_BG   = (255, 249, 230)   # gold tint — lead compound
HEADER_BG = (30,  50,  90)    # deep navy header
GOLD      = (180, 130,   0)   # dark gold (readable on white)
LEAD_BDR  = (210, 155,  10)   # lead border
BLUE_BDR  = (100, 140, 200)   # normal border
TEXT1     = (20,  30,  55)    # near-black
TEXT2     = (60,  80, 120)    # medium navy
TEXT3     = (120, 140, 175)   # muted
GREEN     = (20,  150,  70)   # dark green (readable)
ORANGE    = (200, 100,   0)   # dark orange
DIVIDER   = (220, 225, 235)

def scaffold_info(zid):
    if zid in ('351016','351017','351018','351019'):
        return CARD_BG_A, (160, 100, 0),   "Scaffold A · Isoindolinone-dihydrouracil"
    if zid in ('5024943','5024944','5024945','13662104'):
        return CARD_BG_B, (0, 130, 70),    "Scaffold B · Bicyclic hydrazide-dihydrouracil"
    return     CARD_BG_C, (30, 90, 180),   "Scaffold C · Benzimidazoquinoxaline"

def score_color(s):
    if s <= -8.0: return GREEN
    if s <= -7.5: return ORANGE
    return TEXT2

total_w = PAD + COLS * (MOL_W + PAD)
total_h = HEADER_H + ROWS * (MOL_H + LABEL_H + PAD) + PAD + FOOTER_H

canvas = Image.new('RGB', (total_w, total_h), BG)
draw   = ImageDraw.Draw(canvas)

# Header
draw.rectangle([0, 0, total_w, HEADER_H], fill=HEADER_BG)
draw.line([(0, HEADER_H), (total_w, HEADER_H)], fill=(200, 165, 50), width=3)

# Accent stripe in header
draw.rectangle([0, 0, 6, HEADER_H], fill=(200, 165, 50))

try:
    fnt_title  = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 24)
    fnt_sub    = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",      14)
    fnt_tiny   = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",      12)
    fnt_id     = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 14)
    fnt_score  = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",      13)
    fnt_rank   = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 13)
except:
    fnt_title = fnt_sub = fnt_tiny = fnt_id = fnt_score = fnt_rank = ImageFont.load_default()

draw.text((18, 15), "Pan-Coronavirus RTC Inhibitor Discovery", font=fnt_title, fill=(255, 255, 255))
draw.text((18, 48), "11 Triple-Target Hits  ·  NSP12–NSP7  ·  NSP9–NSP12  ·  NSP12–NSP8  ·  ZINC20 Virtual Screen (9,800 cpds/target)", font=fnt_sub, fill=(180, 200, 235))
draw.text((18, 72), "AutoDock Vina  ·  Exhaustiveness = 16  ·  GIGA-VIN Lab, University of Liège  ·  2026", font=fnt_tiny, fill=(130, 155, 200))

# Badge: 11 hits
badge_x = total_w - 180
draw.rounded_rectangle([badge_x, 20, badge_x+155, 82], radius=8,
                        fill=(50, 75, 130), outline=(200, 165, 50), width=2)
draw.text((badge_x+18, 22), "11 TRIPLE-TARGET", font=fnt_tiny,  fill=(200, 165, 50))
draw.text((badge_x+30, 40), "HITS FOUND",       font=fnt_id,    fill=(255, 255, 255))
draw.text((badge_x+20, 62), "Pan-CoV candidates", font=fnt_tiny, fill=(160, 185, 225))

def render_mol(smi, w, h):
    mol = Chem.MolFromSmiles(smi)
    AllChem.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DCairo(w, h)
    opts = drawer.drawOptions()
    opts.bondLineWidth = 2.0
    opts.padding = 0.12
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return Image.open(io.BytesIO(drawer.GetDrawingText())).convert('RGB')

for idx, (zid, smi, scores, rank) in enumerate(data):
    col = idx % COLS
    row = idx // COLS
    cx  = PAD + col * (MOL_W + PAD)
    cy  = HEADER_H + PAD + row * (MOL_H + LABEL_H + PAD)

    is_lead           = (rank == '#1 Lead')
    sc_bg, sc_col, sc_tag = scaffold_info(zid)
    card_bg           = LEAD_BG if is_lead else sc_bg
    bdr               = LEAD_BDR if is_lead else BLUE_BDR
    bdr_w             = 3 if is_lead else 1

    # Subtle drop shadow
    draw.rounded_rectangle(
        [cx+3, cy+3, cx+MOL_W+3, cy+MOL_H+LABEL_H+3],
        radius=10, fill=(210, 215, 225)
    )
    # Card
    draw.rounded_rectangle(
        [cx, cy, cx+MOL_W, cy+MOL_H+LABEL_H],
        radius=10, fill=card_bg, outline=bdr, width=bdr_w
    )

    # Gold star badge for lead
    if is_lead:
        draw.rounded_rectangle([cx+8, cy+8, cx+70, cy+26], radius=5,
                                fill=(210, 155, 10), outline=None)
        draw.text((cx+14, cy+9), "★ LEAD", font=fnt_rank, fill=(255, 255, 255))

    # Molecule image (white bg — matches light card)
    mol_img = render_mol(smi, MOL_W, MOL_H)
    mask = Image.new('L', (MOL_W, MOL_H), 0)
    mdraw = ImageDraw.Draw(mask)
    mdraw.rounded_rectangle([0, 0, MOL_W, MOL_H + 20], radius=10, fill=255)
    canvas.paste(mol_img, (cx, cy), mask)

    # Divider
    ly = cy + MOL_H
    draw.line([(cx+10, ly), (cx+MOL_W-10, ly)], fill=DIVIDER, width=1)

    # ZINC ID
    label_y = ly + 7
    draw.text((cx+10, label_y), f"ZINC{zid}", font=fnt_id, fill=GOLD if is_lead else TEXT1)
    draw.text((cx+MOL_W-42, label_y), rank, font=fnt_rank, fill=TEXT3)

    # Scores row 1
    sy1 = label_y + 22
    draw.text((cx+10,  sy1), "NSP12-NSP7", font=fnt_score, fill=TEXT3)
    draw.text((cx+108, sy1), f"{scores[0]:.3f}", font=fnt_score, fill=score_color(scores[0]))
    draw.text((cx+185, sy1), "NSP9-NSP12", font=fnt_score, fill=TEXT3)
    draw.text((cx+285, sy1), f"{scores[1]:.3f}", font=fnt_score, fill=score_color(scores[1]))

    # Scores row 2
    sy2 = label_y + 43
    draw.text((cx+10,  sy2), "NSP12-NSP8", font=fnt_score, fill=TEXT3)
    draw.text((cx+108, sy2), f"{scores[2]:.3f}", font=fnt_score, fill=score_color(scores[2]))
    sum_s = sum(scores)
    draw.text((cx+185, sy2), "Σ  score", font=fnt_score, fill=TEXT3)
    sum_col = GOLD if sum_s <= -23.0 else TEXT2
    draw.text((cx+260, sy2), f"{sum_s:.3f}", font=fnt_score, fill=sum_col)

    # Scaffold tag
    sy3 = label_y + 66
    draw.rounded_rectangle([cx+8, sy3-2, cx+MOL_W-8, sy3+16], radius=4, fill=DIVIDER)
    draw.text((cx+14, sy3), sc_tag, font=fnt_tiny, fill=sc_col)

# Footer
fy = total_h - FOOTER_H + 10
draw.line([(0, fy-6), (total_w, fy-6)], fill=DIVIDER, width=1)
draw.text((PAD, fy), "Score key:", font=fnt_tiny, fill=TEXT3)
for x, col, label in [
    (PAD+80,  GREEN,  "≤ −8.0  strong hit"),
    (PAD+215, ORANGE, "≤ −7.5  moderate hit"),
    (PAD+375, TEXT2,  "≤ −7.0  hit threshold"),
]:
    draw.rectangle([x, fy+3, x+13, fy+15], fill=col)
    draw.text((x+17, fy), label, font=fnt_tiny, fill=TEXT1)

draw.text((PAD+565, fy), "★ Lead: ZINC351017  ·  Σ = −24.044 kcal/mol  ·  Scaffold A", font=fnt_tiny, fill=GOLD)

# Scaffold legend
fy2 = fy + 22
draw.text((PAD, fy2), "Scaffolds:", font=fnt_tiny, fill=TEXT3)
for x, bg, tc, label in [
    (PAD+75,  CARD_BG_A, (160,100,0),  "A · Isoindolinone (351016–351019)"),
    (PAD+310, CARD_BG_B, (0,130,70),   "B · Bicyclic hydrazide (5024943–5024945, 13662104)"),
    (PAD+625, CARD_BG_C, (30,90,180),  "C · Benzimidazoquinoxaline (13633805, 13633807, 3169307)"),
]:
    draw.rounded_rectangle([x, fy2, x+13, fy2+14], radius=3, fill=bg, outline=tc, width=1)
    draw.text((x+17, fy2), label, font=fnt_tiny, fill=tc)

out = 'figures/script16/triple_target_hits_structures.png'
canvas.save(out, dpi=(150, 150))
print(f"Saved: {out}")

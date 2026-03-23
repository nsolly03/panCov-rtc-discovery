#!/usr/bin/env python3
"""
CGCP Phase 2 Step 2 - Contact Mapping: NSP12-NSP8
Maps interface residues and computes composite scores.
Anchor: LYS332 (NSP12, cons=1.000)

Run:
    cd ~/projects/rtc-pan-coronavirus
    python3 CGCP/scripts/phase-2/step02_contact_mapping_NSP12-NSP8.py
"""

import os, json, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings('ignore')
from Bio import PDB

plt.rcParams.update({
    'font.family':        'sans-serif',
    'font.sans-serif':    ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size':          9,
    'axes.linewidth':     0.75,
    'axes.spines.top':    False,
    'axes.spines.right':  False,
    'axes.facecolor':     'white',
    'axes.labelsize':     9,
    'axes.labelpad':      6,
    'axes.titlepad':      8,
    'xtick.direction':    'out',
    'ytick.direction':    'out',
    'xtick.major.size':   4,
    'ytick.major.size':   4,
    'xtick.major.width':  0.75,
    'ytick.major.width':  0.75,
    'axes.grid':          False,
    'figure.facecolor':   'white',
    'savefig.dpi':        300,
    'savefig.bbox':       'tight',
    'savefig.facecolor':  'white',
    'savefig.pad_inches': 0.15,
    'pdf.fonttype':       42,
})

BASE     = os.path.expanduser('~/projects/rtc-pan-coronavirus')
PDB_FILE = os.path.join(BASE,
    '03-virtual-screening/NSP12-NSP8_4/receptor_NSP12-NSP8_4.pdb')
OUT_DIR  = os.path.join(BASE,
    'CGCP/02-deep-dive/NSP12-NSP8/step-02-contacts')
os.makedirs(OUT_DIR, exist_ok=True)

CONTACT_CUTOFF = 5.0   # Angstrom
ANCHOR_RES     = 332
ANCHOR_CHAIN   = 'A'

# Real per-residue conservation from aligned FASTA (5 coronaviruses)
NSP12_CONSERVATION = {
    1: 1.0,
    2: 0.8,
    3: 0.6,
    4: 1.0,
    5: 1.0,
    6: 0.4,
    7: 0.6,
    8: 0.4,
    9: 1.0,
    10: 0.4,
    11: 0.4,
    12: 1.0,
    13: 0.6,
    14: 0.8,
    15: 0.8,
    16: 0.8,
    17: 0.8,
    18: 1.0,
    19: 1.0,
    20: 0.6,
    21: 0.8,
    22: 0.4,
    23: 0.8,
    24: 0.6,
    25: 0.6,
    26: 0.4,
    27: 0.8,
    28: 1.0,
    29: 0.4,
    30: 0.4,
    31: 0.4,
    32: 1.0,
    33: 0.8,
    34: 1.0,
    35: 0.8,
    36: 0.6,
    37: 0.6,
    38: 0.8,
    39: 0.4,
    40: 1.0,
    41: 0.8,
    42: 1.0,
    43: 0.4,
    44: 1.0,
    45: 1.0,
    46: 1.0,
    47: 1.0,
    48: 0.4,
    49: 0.6,
    50: 0.4,
    51: 0.4,
    52: 0.6,
    53: 0.4,
    54: 0.4,
    55: 0.4,
    56: 0.4,
    57: 0.6,
    58: 0.6,
    59: 0.6,
    60: 0.8,
    61: 1.0,
    62: 0.8,
    63: 0.6,
    64: 1.0,
    65: 1.0,
    66: 0.8,
    67: 1.0,
    68: 1.0,
    69: 1.0,
    70: 0.6,
    71: 1.0,
    72: 1.0,
    73: 0.8,
    74: 0.8,
    75: 1.0,
    76: 0.8,
    77: 1.0,
    78: 0.6,
    79: 1.0,
    80: 1.0,
    81: 0.8,
    82: 0.5,
    83: 0.8,
    84: 0.8,
    85: 0.6,
    86: 0.6,
    87: 0.8,
    88: 1.0,
    89: 0.8,
    90: 0.8,
    91: 1.0,
    92: 1.0,
    93: 0.4,
    94: 0.6,
    95: 0.6,
    96: 1.0,
    97: 1.0,
    98: 0.8,
    99: 0.4,
    100: 0.6,
    101: 1.0,
    102: 1.0,
    103: 0.4,
    104: 0.4,
    105: 1.0,
    106: 0.4,
    107: 1.0,
    108: 0.6,
    109: 1.0,
    110: 0.6,
    111: 0.8,
    112: 1.0,
    113: 1.0,
    114: 0.6,
    115: 0.8,
    116: 0.8,
    117: 0.4,
    118: 0.8,
    119: 1.0,
    120: 1.0,
    121: 0.8,
    122: 0.4,
    123: 1.0,
    124: 1.0,
    125: 0.6,
    126: 0.6,
    127: 1.0,
    128: 0.6,
    129: 0.6,
    130: 0.8,
    131: 0.4,
    132: 1.0,
    133: 0.8,
    134: 0.6,
    135: 1.0,
    136: 0.4,
    137: 0.8,
    138: 1.0,
    139: 1.0,
    140: 0.6,
    141: 1.0,
    142: 0.6,
    143: 0.2,
    144: 0.6,
    145: 1.0,
    146: 0.4,
    147: 1.0,
    148: 1.0,
    149: 1.0,
    150: 1.0,
    151: 1.0,
    152: 1.0,
    153: 1.0,
    154: 1.0,
    155: 1.0,
    156: 1.0,
    157: 1.0,
    158: 1.0,
    159: 1.0,
    160: 0.8,
    161: 0.6,
    162: 1.0,
    163: 1.0,
    164: 1.0,
    165: 1.0,
    166: 1.0,
    167: 1.0,
    168: 0.8,
    169: 0.4,
    170: 0.4,
    171: 0.2,
    172: 1.0,
    173: 0.8,
    174: 0.4,
    175: 1.0,
    176: 0.8,
    177: 0.8,
    178: 0.4,
    179: 0.6,
    180: 0.6,
    181: 1.0,
    182: 1.0,
    183: 1.0,
    184: 1.0,
    185: 0.6,
    186: 0.6,
    187: 1.0,
    188: 1.0,
    189: 0.6,
    190: 0.6,
    191: 0.4,
    192: 0.6,
    193: 1.0,
    194: 0.4,
    195: 0.6,
    196: 1.0,
    197: 0.8,
    198: 0.6,
    199: 1.0,
    200: 0.4,
    201: 0.6,
    202: 0.4,
    203: 0.6,
    204: 0.4,
    205: 1.0,
    206: 0.4,
    207: 1.0,
    208: 0.4,
    209: 0.6,
    210: 0.6,
    211: 0.4,
    212: 0.6,
    213: 0.4,
    214: 0.6,
    215: 0.8,
    216: 1.0,
    217: 0.8,
    218: 0.8,
    219: 1.0,
    220: 1.0,
    221: 1.0,
    222: 1.0,
    223: 0.8,
    224: 0.4,
    225: 0.6,
    226: 0.4,
    227: 0.2,
    228: 1.0,
    229: 1.0,
    230: 0.4,
    231: 0.6,
    232: 1.0,
    233: 1.0,
    234: 1.0,
    235: 0.8,
    236: 1.0,
    237: 0.6,
    238: 0.8,
    239: 0.6,
    240: 1.0,
    241: 1.0,
    242: 0.8,
    243: 0.8,
    244: 1.0,
    245: 0.6,
    246: 0.6,
    247: 1.0,
    248: 0.4,
    249: 1.0,
    250: 0.8,
    251: 0.6,
    252: 1.0,
    253: 0.6,
    254: 0.8,
    255: 1.0,
    256: 1.0,
    257: 0.8,
    258: 1.0,
    259: 1.0,
    260: 1.0,
    261: 0.6,
    262: 1.0,
    263: 1.0,
    264: 0.6,
    265: 0.8,
    266: 0.6,
    267: 0.4,
    268: 1.0,
    269: 0.4,
    270: 1.0,
    271: 0.4,
    272: 1.0,
    273: 1.0,
    274: 1.0,
    275: 0.8,
    276: 0.6,
    277: 1.0,
    278: 1.0,
    279: 0.6,
    280: 1.0,
    281: 0.6,
    282: 1.0,
    283: 1.0,
    284: 1.0,
    285: 1.0,
    286: 0.6,
    287: 1.0,
    288: 0.6,
    289: 0.6,
    290: 0.4,
    291: 1.0,
    292: 1.0,
    293: 1.0,
    294: 0.8,
    295: 0.6,
    296: 0.6,
    297: 1.0,
    298: 1.0,
    299: 0.6,
    300: 1.0,
    301: 0.4,
    302: 1.0,
    303: 0.4,
    304: 1.0,
    305: 1.0,
    306: 0.8,
    307: 0.6,
    308: 1.0,
    309: 0.8,
    310: 0.4,
    311: 1.0,
    312: 1.0,
    313: 0.6,
    314: 0.4,
    315: 0.6,
    316: 1.0,
    317: 1.0,
    318: 0.8,
    319: 0.4,
    320: 0.6,
    321: 0.6,
    322: 0.6,
    323: 1.0,
    324: 1.0,
    325: 0.6,
    326: 0.6,
    327: 0.6,
    328: 0.4,
    329: 1.0,
    330: 1.0,
    331: 0.6,
    332: 0.6,
    333: 0.6,
    334: 0.8,
    335: 0.6,
    336: 1.0,
    337: 0.6,
    338: 1.0,
    339: 1.0,
    340: 0.4,
    341: 1.0,
    342: 1.0,
    343: 1.0,
    344: 1.0,
    345: 1.0,
    346: 1.0,
    347: 1.0,
    348: 0.6,
    349: 0.6,
    350: 0.6,
    351: 0.6,
    352: 0.6,
    353: 0.6,
    354: 1.0,
    355: 0.8,
    356: 1.0,
    357: 0.8,
    358: 1.0,
    359: 1.0,
    360: 0.6,
    361: 1.0,
    362: 1.0,
    363: 0.6,
    364: 0.6,
    365: 1.0,
    366: 1.0,
    367: 0.8,
    368: 1.0,
    369: 0.4,
    370: 0.6,
    371: 0.8,
    372: 0.6,
    373: 1.0,
    374: 1.0,
    375: 1.0,
    376: 0.6,
    377: 1.0,
    378: 1.0,
    379: 1.0,
    380: 0.6,
    381: 0.6,
    382: 0.6,
    383: 1.0,
    384: 1.0,
    385: 1.0,
    386: 1.0,
    387: 1.0,
    388: 1.0,
    389: 0.8,
    390: 1.0,
    391: 0.6,
    392: 1.0,
    393: 0.6,
    394: 1.0,
    395: 1.0,
    396: 0.8,
    397: 0.4,
    398: 1.0,
    399: 0.6,
    400: 0.8,
    401: 0.8,
    402: 1.0,
    403: 0.8,
    404: 1.0,
    405: 1.0,
    406: 0.6,
    407: 1.0,
    408: 1.0,
    409: 0.8,
    410: 0.6,
    411: 1.0,
    412: 1.0,
    413: 0.4,
    414: 1.0,
    415: 0.4,
    416: 0.6,
    417: 0.6,
    418: 0.4,
    419: 0.4,
    420: 0.6,
    421: 0.8,
    422: 0.8,
    423: 0.4,
    424: 0.6,
    425: 1.0,
    426: 1.0,
    427: 0.8,
    428: 0.6,
    429: 1.0,
    430: 0.6,
    431: 1.0,
    432: 1.0,
    433: 1.0,
    434: 0.8,
    435: 0.6,
    436: 0.8,
    437: 0.4,
    438: 0.6,
    439: 1.0,
    440: 0.6,
    441: 1.0,
    442: 0.6,
    443: 1.0,
    444: 1.0,
    445: 0.6,
    446: 1.0,
    447: 1.0,
    448: 1.0,
    449: 1.0,
    450: 0.4,
    451: 1.0,
    452: 0.6,
    453: 1.0,
    454: 1.0,
    455: 0.6,
    456: 1.0,
    457: 1.0,
    458: 1.0,
    459: 0.6,
    460: 0.8,
    461: 1.0,
    462: 1.0,
    463: 0.6,
    464: 1.0,
    465: 0.6,
    466: 1.0,
    467: 1.0,
    468: 0.8,
    469: 0.6,
    470: 1.0,
    471: 1.0,
    472: 0.8,
    473: 0.8,
    474: 1.0,
    475: 0.8,
    476: 0.4,
    477: 1.0,
    478: 1.0,
    479: 1.0,
    480: 1.0,
    481: 0.8,
    482: 0.6,
    483: 1.0,
    484: 1.0,
    485: 0.6,
    486: 1.0,
    487: 1.0,
    488: 0.6,
    489: 1.0,
    490: 1.0,
    491: 1.0,
    492: 1.0,
    493: 1.0,
    494: 1.0,
    495: 1.0,
    496: 0.6,
    497: 1.0,
    498: 0.6,
    499: 1.0,
    500: 1.0,
    501: 1.0,
    502: 1.0,
    503: 1.0,
    504: 0.6,
    505: 1.0,
    506: 1.0,
    507: 1.0,
    508: 0.6,
    509: 0.6,
    510: 0.8,
    511: 1.0,
    512: 1.0,
    513: 1.0,
    514: 0.6,
    515: 1.0,
    516: 1.0,
    517: 0.6,
    518: 1.0,
    519: 1.0,
    520: 1.0,
    521: 0.4,
    522: 1.0,
    523: 1.0,
    524: 1.0,
    525: 0.8,
    526: 0.6,
    527: 0.8,
    528: 1.0,
    529: 1.0,
    530: 0.6,
    531: 1.0,
    532: 1.0,
    533: 0.8,
    534: 1.0,
    535: 1.0,
    536: 1.0,
    537: 1.0,
    538: 0.6,
    539: 1.0,
    540: 1.0,
    541: 1.0,
    542: 1.0,
    543: 1.0,
    544: 1.0,
    545: 0.4,
    546: 0.8,
    547: 1.0,
    548: 1.0,
    549: 0.8,
    550: 0.8,
    551: 0.6,
    552: 0.6,
    553: 0.4,
    554: 0.8,
    555: 1.0,
    556: 0.6,
    557: 0.6,
    558: 1.0,
    559: 0.6,
    560: 1.0,
    561: 1.0,
    562: 1.0,
    563: 1.0,
    564: 1.0,
    565: 1.0,
    566: 1.0,
    567: 1.0,
    568: 1.0,
    569: 1.0,
    570: 1.0,
    571: 1.0,
    572: 0.8,
    573: 1.0,
    574: 0.8,
    575: 1.0,
    576: 0.4,
    577: 1.0,
    578: 0.6,
    579: 0.4,
    580: 0.6,
    581: 0.6,
    582: 0.6,
    583: 0.6,
    584: 1.0,
    585: 0.6,
    586: 0.6,
    587: 1.0,
    588: 1.0,
    589: 0.4,
    590: 0.8,
    591: 1.0,
    592: 1.0,
    593: 0.6,
    594: 0.4,
    595: 0.6,
    596: 0.6,
    597: 0.8,
    598: 1.0,
    599: 1.0,
    600: 1.0,
    601: 1.0,
    602: 0.6,
    603: 1.0,
    604: 1.0,
    605: 0.6,
    606: 1.0,
    607: 1.0,
    608: 1.0,
    609: 1.0,
    610: 0.6,
    611: 1.0,
    612: 0.4,
    613: 1.0,
    614: 0.4,
    615: 0.6,
    616: 0.6,
    617: 1.0,
    618: 0.6,
    619: 0.4,
    620: 1.0,
    621: 0.6,
    622: 1.0,
    623: 1.0,
    624: 1.0,
    625: 1.0,
    626: 1.0,
    627: 0.6,
    628: 1.0,
    629: 1.0,
    630: 1.0,
    631: 1.0,
    632: 0.8,
    633: 1.0,
    634: 1.0,
    635: 1.0,
    636: 1.0,
    637: 1.0,
    638: 1.0,
    639: 0.8,
    640: 1.0,
    641: 1.0,
    642: 1.0,
    643: 0.4,
    644: 1.0,
    645: 1.0,
    646: 0.8,
    647: 0.6,
    648: 0.6,
    649: 1.0,
    650: 0.6,
    651: 0.8,
    652: 0.6,
    653: 0.8,
    654: 0.8,
    655: 0.8,
    656: 0.4,
    657: 0.4,
    658: 0.6,
    659: 0.6,
    660: 0.6,
    661: 0.6,
    662: 0.4,
    663: 0.6,
    664: 0.6,
    665: 0.4,
    666: 1.0,
    667: 0.6,
    668: 0.4,
    669: 0.8,
    670: 1.0,
    671: 0.4,
    672: 0.6,
    673: 1.0,
    674: 1.0,
    675: 0.4,
    676: 0.6,
    677: 0.4,
    678: 1.0,
    679: 1.0,
    680: 0.6,
    681: 0.4,
    682: 0.4,
    683: 0.8,
    684: 0.8,
    685: 0.4,
    686: 0.4,
    687: 1.0,
    688: 0.8,
    689: 0.8,
    690: 0.4,
    691: 0.6,
    692: 1.0,
    693: 0.6,
    694: 0.8,
    695: 1.0,
    696: 0.6,
    697: 1.0,
    698: 1.0,
    699: 1.0,
    700: 1.0,
    701: 1.0,
    702: 1.0,
    703: 1.0,
    704: 1.0,
    705: 1.0,
    706: 1.0,
    707: 1.0,
    708: 0.4,
    709: 1.0,
    710: 1.0,
    711: 1.0,
    712: 0.8,
    713: 1.0,
    714: 0.6,
    715: 0.4,
    716: 1.0,
    717: 1.0,
    718: 0.4,
    719: 0.4,
    720: 1.0,
    721: 0.6,
    722: 0.6,
    723: 1.0,
    724: 0.4,
    725: 1.0,
    726: 0.4,
    727: 0.6,
    728: 1.0,
    729: 1.0,
    730: 0.6,
    731: 0.6,
    732: 1.0,
    733: 1.0,
    734: 1.0,
    735: 1.0,
    736: 1.0,
    737: 0.8,
    738: 1.0,
    739: 1.0,
    740: 1.0,
    741: 1.0,
    742: 0.6,
    743: 0.8,
    744: 1.0,
    745: 1.0,
    746: 1.0,
    747: 0.6,
    748: 1.0,
    749: 0.6,
    750: 1.0,
    751: 1.0,
    752: 0.6,
    753: 0.8,
    754: 1.0,
    755: 1.0,
    756: 1.0,
    757: 1.0,
    758: 1.0,
    759: 1.0,
    760: 1.0,
    761: 1.0,
    762: 1.0,
    763: 1.0,
    764: 0.8,
    765: 0.4,
    766: 0.6,
    767: 0.6,
    768: 0.6,
    769: 0.6,
    770: 0.8,
    771: 0.6,
    772: 0.4,
    773: 0.6,
    774: 0.8,
    775: 1.0,
    776: 1.0,
    777: 1.0,
    778: 1.0,
    779: 1.0,
    780: 1.0,
    781: 1.0,
    782: 1.0,
    783: 1.0,
    784: 0.8,
    785: 0.6,
    786: 1.0,
    787: 1.0,
    788: 0.6,
    789: 1.0,
    790: 1.0,
    791: 1.0,
    792: 1.0,
    793: 0.8,
    794: 0.8,
    795: 1.0,
    796: 1.0,
    797: 1.0,
    798: 0.6,
    799: 0.6,
    800: 0.6,
    801: 0.6,
    802: 0.4,
    803: 1.0,
    804: 1.0,
    805: 0.6,
    806: 1.0,
    807: 1.0,
    808: 1.0,
    809: 1.0,
    810: 1.0,
    811: 1.0,
    812: 1.0,
    813: 1.0,
    814: 1.0,
    815: 1.0,
    816: 0.6,
    817: 1.0,
    818: 1.0,
    819: 0.8,
    820: 0.6,
    821: 0.4,
    822: 1.0,
    823: 1.0,
    824: 0.4,
    825: 0.4,
    826: 1.0,
    827: 1.0,
    828: 0.4,
    829: 0.4,
    830: 0.6,
    831: 1.0,
    832: 0.6,
    833: 0.6,
    834: 0.6,
    835: 0.4,
    836: 0.6,
    837: 1.0,
    838: 0.4,
    839: 0.6,
    840: 0.4,
    841: 1.0,
    842: 0.6,
    843: 0.6,
    844: 0.6,
    845: 0.6,
    846: 1.0,
    847: 0.6,
    848: 0.6,
    849: 0.6,
    850: 1.0,
    851: 1.0,
    852: 0.6,
    853: 1.0,
    854: 0.4,
    855: 0.4,
    856: 0.6,
    857: 0.6,
    858: 0.4,
    859: 0.6,
    860: 0.6,
    861: 0.6,
    862: 1.0,
    863: 0.6,
    864: 0.6,
    865: 0.4,
    866: 1.0,
    867: 1.0,
    868: 0.4,
    869: 0.4,
    870: 0.8,
    871: 1.0,
    872: 0.4,
    873: 0.5,
    874: 0.5,
    875: 1.0,
    876: 0.5,
    877: 1.0,
    878: 1.0,
    879: 1.0,
    880: 0.667,
    881: 1.0,
    882: 0.333,
    883: 1.0,
    884: 1.0,
    885: 1.0,
    886: 1.0,
    887: 1.0,
    888: 1.0,
    889: 1.0,
    890: 1.0,
    891: 1.0,
    892: 1.0,
    893: 1.0,
    894: 1.0,
    895: 1.0,
    896: 1.0,
    897: 1.0,
    898: 1.0,
    899: 1.0,
    900: 1.0,
    901: 1.0,
    902: 1.0,
    903: 1.0,
    904: 1.0,
    905: 1.0,
    906: 1.0,
    907: 1.0,
    908: 1.0,
    909: 1.0,
    910: 1.0,
    911: 1.0,
    912: 1.0,
    913: 1.0,
    914: 1.0,
    915: 1.0,
    916: 1.0,
    917: 1.0,
    918: 1.0,
    919: 1.0,
    920: 1.0,
    921: 1.0,
    922: 1.0,
}
NSP8_CONSERVATION = {
    1: 1.0,
    2: 1.0,
    3: 1.0,
    4: 0.6,
    5: 1.0,
    6: 1.0,
    7: 0.6,
    8: 0.6,
    9: 0.6,
    10: 0.6,
    11: 0.6,
    12: 0.6,
    13: 1.0,
    14: 0.8,
    15: 0.4,
    16: 1.0,
    17: 1.0,
    18: 1.0,
    19: 0.6,
    20: 0.6,
    21: 0.6,
    22: 1.0,
    23: 1.0,
    24: 0.6,
    25: 1.0,
    26: 0.8,
    27: 0.8,
    28: 0.6,
    29: 1.0,
    30: 1.0,
    31: 0.6,
    32: 0.6,
    33: 0.4,
    34: 0.6,
    35: 0.6,
    36: 0.4,
    37: 0.8,
    38: 0.6,
    39: 0.8,
    40: 0.6,
    41: 1.0,
    42: 1.0,
    43: 0.6,
    44: 1.0,
    45: 0.6,
    46: 0.8,
    47: 1.0,
    48: 1.0,
    49: 0.6,
    50: 0.6,
    51: 1.0,
    52: 1.0,
    53: 0.6,
    54: 0.8,
    55: 0.4,
    56: 0.8,
    57: 0.6,
    58: 0.8,
    59: 1.0,
    60: 0.4,
    61: 0.8,
    62: 1.0,
    63: 0.4,
    64: 1.0,
    65: 1.0,
    66: 0.4,
    67: 0.2,
    68: 0.4,
    69: 0.6,
    70: 0.4,
    71: 0.4,
    72: 0.4,
    73: 0.4,
    74: 0.4,
    75: 0.4,
    76: 0.4,
    77: 0.4,
    78: 0.4,
    79: 0.4,
    80: 1.0,
    81: 1.0,
    82: 0.6,
    83: 0.6,
    84: 0.6,
    85: 1.0,
    86: 0.4,
    87: 0.6,
    88: 0.6,
    89: 0.4,
    90: 0.8,
    91: 0.4,
    92: 0.8,
    93: 0.4,
    94: 1.0,
    95: 0.4,
    96: 0.6,
    97: 0.6,
    98: 1.0,
    99: 0.4,
    100: 0.6,
    101: 0.6,
    102: 0.6,
    103: 0.6,
    104: 0.4,
    105: 0.5,
    106: 0.6,
    107: 0.8,
    108: 0.4,
    109: 0.6,
    110: 0.6,
    111: 1.0,
    112: 0.4,
    113: 1.0,
    114: 1.0,
    115: 1.0,
    116: 0.8,
    117: 0.6,
    118: 0.6,
    119: 0.6,
    120: 0.4,
    121: 1.0,
    122: 0.6,
    123: 1.0,
    124: 1.0,
    125: 1.0,
    126: 1.0,
    127: 1.0,
    128: 1.0,
    129: 1.0,
    130: 0.75,
    131: 0.667,
    132: 1.0,
    133: 1.0,
    134: 1.0,
    135: 0.5,
    136: 1.0,
    137: 0.5,
    138: 1.0,
    139: 0.5,
    140: 0.5,
    141: 1.0,
    142: 0.5,
    143: 1.0,
    144: 1.0,
    145: 1.0,
    146: 1.0,
    147: 1.0,
    148: 1.0,
    149: 1.0,
    150: 1.0,
    151: 1.0,
    152: 1.0,
    153: 1.0,
    154: 1.0,
    155: 1.0,
    156: 1.0,
    157: 1.0,
    158: 1.0,
    159: 1.0,
    160: 1.0,
    161: 1.0,
    162: 1.0,
    163: 1.0,
    164: 1.0,
    165: 1.0,
    166: 1.0,
    167: 1.0,
    168: 1.0,
    169: 1.0,
    170: 1.0,
    171: 1.0,
    172: 1.0,
    173: 1.0,
    174: 1.0,
    175: 1.0,
    176: 1.0,
    177: 1.0,
    178: 1.0,
    179: 1.0,
    180: 1.0,
    181: 1.0,
    182: 1.0,
    183: 1.0,
    184: 1.0,
    185: 1.0,
    186: 1.0,
    187: 1.0,
    188: 1.0,
    189: 1.0,
    190: 1.0,
    191: 1.0,
    192: 1.0,
    193: 1.0,
}

def get_default_cons(pos, chain='NSP12'):
    if chain == 'NSP12':
        return NSP12_CONSERVATION.get(pos, 0.600)
    return NSP8_CONSERVATION.get(pos, 0.600)



def compute_bsa(res, chain_partner, cutoff=8.0):
    """Approximate BSA as number of partner atoms within cutoff."""
    count = 0
    for atom in res.get_atoms():
        for partner_res in chain_partner.get_residues():
            if not PDB.is_aa(partner_res):
                continue
            for patom in partner_res.get_atoms():
                d = atom - patom
                if d < cutoff:
                    count += 1
    return count


def map_contacts(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('cx', pdb_file)[0]
    chain_a = structure['A']   # NSP12
    chain_b = structure['B']   # NSP8

    # Get all residues
    res_a = [r for r in chain_a.get_residues() if PDB.is_aa(r)]
    res_b = [r for r in chain_b.get_residues() if PDB.is_aa(r)]

    # Find interface residues within CONTACT_CUTOFF
    interface_a = {}
    interface_b = {}

    for ra in res_a:
        pos_a = ra.id[1]
        contacts = 0
        partners = []
        for rb in res_b:
            for atom_a in ra.get_atoms():
                for atom_b in rb.get_atoms():
                    d = atom_a - atom_b
                    if d < CONTACT_CUTOFF:
                        contacts += 1
                        pos_b = rb.id[1]
                        if pos_b not in partners:
                            partners.append(pos_b)
        if contacts > 0:
            ca = ra['CA'].coord if 'CA' in ra else None
            interface_a[pos_a] = {
                'residue':   f'NSP12-{ra.resname}{pos_a}',
                'chain':     'NSP12',
                'position':  pos_a,
                'aa':        ra.resname,
                'contacts':  contacts,
                'partners':  partners,
                'ca_coord':  ca.tolist() if ca is not None else None,
            }

    for rb in res_b:
        pos_b = rb.id[1]
        contacts = 0
        partners = []
        for ra in res_a:
            for atom_b in rb.get_atoms():
                for atom_a in ra.get_atoms():
                    d = atom_b - atom_a
                    if d < CONTACT_CUTOFF:
                        contacts += 1
                        pos_a = ra.id[1]
                        if pos_a not in partners:
                            partners.append(pos_a)
        if contacts > 0:
            ca = rb['CA'].coord if 'CA' in rb else None
            interface_b[pos_b] = {
                'residue':   f'NSP8-{rb.resname}{pos_b}',
                'chain':     'NSP8',
                'position':  pos_b,
                'aa':        rb.resname,
                'contacts':  contacts,
                'partners':  partners,
                'ca_coord':  ca.tolist() if ca is not None else None,
            }

    return interface_a, interface_b


def compute_composite(contacts, bsa_approx, conservation):
    """
    Composite score = weighted combination of:
      - contact score (normalized)
      - BSA approximation (normalized)
      - conservation
    """
    # Will normalize after collecting all values
    return {
        'contact_score': contacts,
        'bsa_approx':    bsa_approx,
        'conservation':  conservation,
    }


def build_records(interface_a, interface_b, chain_a, chain_b):
    records = []

    all_contacts = ([v['contacts'] for v in interface_a.values()] +
                    [v['contacts'] for v in interface_b.values()])
    max_contacts = max(all_contacts) if all_contacts else 1

    for pos, data in interface_a.items():
        cons = get_default_cons(pos)
        norm_contacts = data['contacts'] / max_contacts
        composite = round(
            0.4 * norm_contacts + 0.3 * min(data['contacts']/50, 1.0)
            + 0.3 * cons, 4)
        records.append({
            'residue':       data['residue'],
            'chain':         'NSP12',
            'position':      pos,
            'aa':            data['aa'],
            'contact_score': data['contacts'],
            'partners':      ','.join(map(str, data['partners'])),
            'n_partners':    len(data['partners']),
            'conservation':  cons,
            'composite':     composite,
            'is_anchor':     1 if pos == ANCHOR_RES else 0,
            'ca_x': round(data['ca_coord'][0], 3)
                    if data['ca_coord'] else None,
            'ca_y': round(data['ca_coord'][1], 3)
                    if data['ca_coord'] else None,
            'ca_z': round(data['ca_coord'][2], 3)
                    if data['ca_coord'] else None,
        })

    for pos, data in interface_b.items():
        cons = get_default_cons(pos)
        norm_contacts = data['contacts'] / max_contacts
        composite = round(
            0.4 * norm_contacts + 0.3 * min(data['contacts']/50, 1.0)
            + 0.3 * cons, 4)
        records.append({
            'residue':       data['residue'],
            'chain':         'NSP8',
            'position':      pos,
            'aa':            data['aa'],
            'contact_score': data['contacts'],
            'partners':      ','.join(map(str, data['partners'])),
            'n_partners':    len(data['partners']),
            'conservation':  cons,
            'composite':     composite,
            'is_anchor':     0,
            'ca_x': round(data['ca_coord'][0], 3)
                    if data['ca_coord'] else None,
            'ca_y': round(data['ca_coord'][1], 3)
                    if data['ca_coord'] else None,
            'ca_z': round(data['ca_coord'][2], 3)
                    if data['ca_coord'] else None,
        })

    # Normalize composite to 0-1
    max_comp = max(r['composite'] for r in records) if records else 1
    for r in records:
        r['composite'] = round(r['composite'] / max_comp, 4)

    # Sort: anchor first, then by composite
    records.sort(key=lambda x: (
        0 if x['is_anchor'] else 1,
        -x['composite']
    ))

    return records


def save_outputs(records):
    # TSV
    tsv_path = os.path.join(OUT_DIR,
        'contact_map_NSP12-NSP8.tsv')
    fields = ['residue','chain','position','aa',
              'contact_score','partners','n_partners',
              'conservation','composite','is_anchor',
              'ca_x','ca_y','ca_z']
    with open(tsv_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields,
                           delimiter='\t')
        w.writeheader()
        w.writerows(records)
    print(f"  TSV: {tsv_path}")

    # JSON
    json_path = os.path.join(OUT_DIR,
        'contact_map_NSP12-NSP8.json')
    summary = {
        'interface':     'NSP12-NSP8',
        'pdb':           PDB_FILE,
        'anchor':        'LYS332',
        'cutoff_A':      CONTACT_CUTOFF,
        'n_interface':   len(records),
        'n_nsp12':       sum(1 for r in records
                             if r['chain']=='NSP12'),
        'n_nsp8':        sum(1 for r in records
                             if r['chain']=='NSP8'),
        'top5': [{'residue': r['residue'],
                  'composite': r['composite'],
                  'conservation': r['conservation']}
                 for r in records[:5]],
        'records': records,
    }
    with open(json_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"  JSON: {json_path}")
    return summary


def make_figure(records, summary):
    fig, axes = plt.subplots(1, 2,
                             figsize=(13, 5.5))
    ax1, ax2 = axes

    # Colors
    NSP12_COL = '#D7191C'
    NSP12_FILL = '#FDAE61'
    NSP8_COL  = '#2166AC'
    NSP8_FILL  = '#92C5DE'
    ANCHOR_COL = '#1A1A1A'

    # ── Panel A: Composite score per residue ──────────────────
    top = records[:30]  # top 30 by composite
    labels = [r['aa'][:3]+str(r['position'])
              for r in top]
    vals   = [r['composite'] for r in top]
    cols   = [ANCHOR_COL if r['is_anchor']
              else (NSP12_COL if r['chain']=='NSP12'
                    else NSP8_COL)
              for r in top]
    fills  = [ANCHOR_COL if r['is_anchor']
              else (NSP12_FILL if r['chain']=='NSP12'
                    else NSP8_FILL)
              for r in top]

    x = np.arange(len(top))
    for xi, val, fc, ec in zip(x, vals, fills, cols):
        ax1.bar(xi, val, width=0.42,
                facecolor=fc, edgecolor=ec,
                linewidth=0.9, zorder=3,
                clip_on=False)
    ax1.scatter(x, vals,
                color=cols, s=16, zorder=5,
                clip_on=False)

    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, fontsize=6.5,
                        rotation=90, ha='center',
                        va='top')
    ax1.tick_params(axis='x', pad=2,
                    length=3, width=0.75)
    ax1.set_xlim(-0.6, len(top)-0.4)
    ax1.set_ylabel('CGCP composite score', fontsize=9)
    ax1.spines['left'].set_position(('outward', 5))
    ax1.spines['bottom'].set_position(('outward', 5))
    ax1.spines['left'].set_bounds(0, 1.0)
    ax1.set_ylim(-0.04, 1.2)
    ax1.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])

    handles = [
        mpatches.Patch(facecolor=ANCHOR_COL,
                       edgecolor=ANCHOR_COL,
                       label='LYS332 anchor'),
        mpatches.Patch(facecolor=NSP12_FILL,
                       edgecolor=NSP12_COL,
                       label='NSP12'),
        mpatches.Patch(facecolor=NSP8_FILL,
                       edgecolor=NSP8_COL,
                       label='NSP8'),
    ]
    ax1.legend(handles=handles, fontsize=7.5,
               frameon=True, facecolor='white',
               edgecolor='#CCCCCC')
    ax1.set_title(
        'A   Contact map composite scores (top 30)\n'
        f"NSP12={summary['n_nsp12']} res | "
        f"NSP8={summary['n_nsp8']} res | "
        f"total={summary['n_interface']}",
        loc='left', fontsize=8.5, pad=6,
        linespacing=1.5)

    # ── Panel B: Conservation vs contact score ────────────────
    for r in records:
        col  = (ANCHOR_COL if r['is_anchor']
                else (NSP12_COL if r['chain']=='NSP12'
                      else NSP8_COL))
        fill = (ANCHOR_COL if r['is_anchor']
                else (NSP12_FILL if r['chain']=='NSP12'
                      else NSP8_FILL))
        size = 120 if r['is_anchor'] else 35
        marker = '*' if r['is_anchor'] else 'o'
        ax2.scatter(r['contact_score'],
                    r['conservation'],
                    color=fill, edgecolors=col,
                    linewidths=0.7, s=size,
                    marker=marker, alpha=0.8,
                    zorder=4)
        if r['is_anchor'] or r['composite'] > 0.7:
            ax2.annotate(
                r['aa'][:3]+str(r['position']),
                (r['contact_score'],
                 r['conservation']),
                fontsize=6.5, color=col,
                xytext=(4, 2),
                textcoords='offset points')

    ax2.set_xlabel('Contact score (atom pairs)',
                   fontsize=9)
    ax2.set_ylabel('Conservation score', fontsize=9)
    ax2.spines['left'].set_position(('outward', 5))
    ax2.spines['bottom'].set_position(('outward', 5))
    ax2.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax2.set_ylim(-0.05, 1.1)
    ax2.set_title(
        'B   Conservation vs contact score\n'
        '(star = LYS332 anchor)',
        loc='left', fontsize=8.5, pad=6,
        linespacing=1.5)

    fig.tight_layout(pad=2.0)
    path = os.path.join(OUT_DIR,
        'Fig_Step02_ContactMap_NSP12-NSP8.png')
    fig.savefig(path)
    plt.close()
    print(f"  Figure: {path}")


def print_summary(records, summary):
    print('\n' + '='*60)
    print('STEP 2 CONTACT MAPPING — NSP12-NSP8')
    print('='*60)
    print(f"  Total interface residues: {summary['n_interface']}")
    print(f"  NSP12: {summary['n_nsp12']} | "
          f"NSP8: {summary['n_nsp8']}")
    print(f"\n  Top 10 by composite score:")
    print(f"  {'Residue':<20} {'Chain':<7} "
          f"{'Contacts':>8} {'Cons':>6} {'Comp':>6}")
    print(f"  {'-'*55}")
    for r in records[:10]:
        sym = '★' if r['is_anchor'] else ' '
        print(f"  {sym} {r['aa']}{r['position']:<15} "
              f"{r['chain']:<7} "
              f"{r['contact_score']:>8} "
              f"{r['conservation']:>6.3f} "
              f"{r['composite']:>6.3f}")
    print('='*60)


if __name__ == '__main__':
    print('CGCP Phase 2 Step 2 — NSP12-NSP8 Contact Mapping')

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('cx', PDB_FILE)[0]
    chain_a = structure['A']
    chain_b = structure['B']

    print('  Mapping interface contacts...')
    interface_a, interface_b = map_contacts(PDB_FILE)
    records = build_records(interface_a, interface_b,
                            chain_a, chain_b)

    summary = save_outputs(records)
    make_figure(records, summary)
    print_summary(records, summary)

    print(f"\nOutputs: {OUT_DIR}")
    print("Next: Step 3 — feature classification")

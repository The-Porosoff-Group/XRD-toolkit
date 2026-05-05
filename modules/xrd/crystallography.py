"""
modules/xrd/crystallography.py
Pure-Python crystallography engine.
No pymatgen dependency — implements everything needed for Le Bail refinement.
"""

import math, re, itertools
import numpy as np


# ─────────────────────────────────────────────────────────────────────────────
# UNIT CELL
# ─────────────────────────────────────────────────────────────────────────────

CRYSTAL_SYSTEMS = {
    'cubic':        lambda a,b,c,al,be,ga: (a, a, a, 90, 90, 90),
    'tetragonal':   lambda a,b,c,al,be,ga: (a, a, c, 90, 90, 90),
    'orthorhombic': lambda a,b,c,al,be,ga: (a, b, c, 90, 90, 90),
    'hexagonal':    lambda a,b,c,al,be,ga: (a, a, c, 90, 90,120),
    'trigonal':     lambda a,b,c,al,be,ga: (a, a, c, 90, 90,120),
    'monoclinic':   lambda a,b,c,al,be,ga: (a, b, c, 90, be, 90),
    'triclinic':    lambda a,b,c,al,be,ga: (a, b, c, al, be, ga),
}

def cell_volume(a, b, c, al, be, ga):
    """Unit cell volume in Å³."""
    al_r = math.radians(al)
    be_r = math.radians(be)
    ga_r = math.radians(ga)
    return a*b*c*math.sqrt(
        1 - math.cos(al_r)**2 - math.cos(be_r)**2 - math.cos(ga_r)**2
        + 2*math.cos(al_r)*math.cos(be_r)*math.cos(ga_r)
    )

def d_spacing(h, k, l, a, b, c, al, be, ga, system):
    """Calculate d-spacing for given hkl and cell parameters."""
    # Guard against None cell parameters
    if any(v is None for v in [a, b, c, al, be, ga]):
        return None
    a = float(a); b = float(b); c = float(c)
    al = float(al); be = float(be); ga = float(ga)
    system = system.lower()
    try:
        if system in ('cubic',):
            return a / math.sqrt(h**2 + k**2 + l**2)
        elif system in ('tetragonal',):
            return 1 / math.sqrt((h**2+k**2)/a**2 + l**2/c**2)
        elif system in ('orthorhombic',):
            return 1 / math.sqrt(h**2/a**2 + k**2/b**2 + l**2/c**2)
        elif system in ('hexagonal', 'trigonal'):
            return 1 / math.sqrt((4/3)*(h**2+h*k+k**2)/a**2 + l**2/c**2)
        elif system in ('monoclinic',):
            be_r = math.radians(be)
            sin_be = math.sin(be_r)
            cos_be = math.cos(be_r)
            return 1 / math.sqrt(
                (1/sin_be**2) * (h**2/a**2 + k**2*sin_be**2/b**2 + l**2/c**2
                                 - 2*h*l*cos_be/(a*c))
            )
        elif system in ('triclinic',):
            al_r, be_r, ga_r = math.radians(al), math.radians(be), math.radians(ga)
            V = cell_volume(a, b, c, al, be, ga)
            S11 = b**2*c**2*math.sin(al_r)**2
            S22 = a**2*c**2*math.sin(be_r)**2
            S33 = a**2*b**2*math.sin(ga_r)**2
            S12 = a*b*c**2*(math.cos(al_r)*math.cos(be_r)-math.cos(ga_r))
            S23 = a**2*b*c*(math.cos(be_r)*math.cos(ga_r)-math.cos(al_r))
            S13 = a*b**2*c*(math.cos(ga_r)*math.cos(al_r)-math.cos(be_r))
            return V / math.sqrt(S11*h**2 + S22*k**2 + S33*l**2
                                  + 2*S12*h*k + 2*S23*k*l + 2*S13*h*l)
        else:
            # Fallback: orthorhombic
            return 1 / math.sqrt(h**2/a**2 + k**2/b**2 + l**2/c**2)
    except (ValueError, ZeroDivisionError):
        return None


# ─────────────────────────────────────────────────────────────────────────────
# SYSTEMATIC ABSENCE RULES
# ─────────────────────────────────────────────────────────────────────────────

def is_allowed(h, k, l, spacegroup_number):
    """
    Systematic absence filter based on lattice centering.
    Covers the most common catalyst phases correctly.
    """
    sg = spacegroup_number

    # Always remove (000)
    if h == 0 and k == 0 and l == 0:
        return False

    # ── Lattice centering conditions ──────────────────────────────────────────

    # F-centred (FCC): h,k,l must be all odd or all even
    # F-lattice cubic SGs: 196,202,203,209,210,216,219,225,226,227,228
    _F_cubic = {196,202,203,209,210,216,219,225,226,227,228}
    if sg in _F_cubic:
        if not ((h%2 == k%2 == l%2 == 0) or (h%2==1 and k%2==1 and l%2==1)):
            return False

    # I-centred (BCC): h+k+l must be even
    # I-lattice cubic SGs: 197,199,204,206,211,214,217,220,229,230
    _I_cubic = {197,199,204,206,211,214,217,220,229,230}
    if sg in _I_cubic:
        if (h+k+l) % 2 != 0:
            return False

    # I-centred tetragonal: h+k+l even (SGs 79-88, 97-98, 107-110, 119-122, 139-142)
    _I_tetrag = set(range(79,89)) | set(range(97,99)) | set(range(107,111)) | \
                set(range(119,123)) | set(range(139,143))
    if sg in _I_tetrag:
        if (h+k+l) % 2 != 0:
            return False

    # I-centred orthorhombic: h+k+l even (SGs 23-24, 44-46, 71-74)
    _I_ortho = {23,24,44,45,46,71,72,73,74}
    if sg in _I_ortho:
        if (h+k+l) % 2 != 0:
            return False

    # C-centred: h+k even
    # NB: Only include actually C-centred SGs. SG 13 (P2/c) and 14 (P21/c)
    # are P-lattice despite being in the monoclinic range 12–15.
    # SGs 38–41 are A-centred (k+l=2n), not C-centred.
    _C_centred = ({5} | {8, 9} | {12, 15} |
                  {20, 21} | {35, 36, 37} | set(range(63, 69)))
    if sg in _C_centred:
        if (h+k) % 2 != 0:
            return False

    # A-centred: k+l even (SGs 38–41: Amm2, Aem2, Ama2, Aea2)
    _A_centred = {38, 39, 40, 41}
    if sg in _A_centred:
        if (k+l) % 2 != 0:
            return False

    # ── Screw axis / glide plane conditions ────────────────────────────────────
    #
    # For each space group we enforce the ITC Vol A "reflection conditions"
    # (equivalent of "systematic absences").  Missing rules here cause
    # spurious zero-intensity sticks in the diffraction pattern.
    #
    # The rules below cover: (a) all cubic SGs already handled above via
    # lattice centering, plus specific glide rules; (b) common monoclinic,
    # orthorhombic, tetragonal, and hexagonal SGs encountered in catalysis.

    # ── Monoclinic (non-duplicated SGs) ──────────────────────────────────

    # P21 (#4): 0k0: k = 2n
    if sg == 4:
        if h == 0 and l == 0 and k % 2 != 0: return False

    # ── Orthorhombic (non-duplicated SGs) ─────────────────────────────────

    # Pna21 (#33): n-glide ⊥ a, a-glide ⊥ b, 21 screw ∥ c
    #   0kl: k+l = 2n (n-glide)
    #   h0l: h = 2n   (a-glide)
    #   00l: l = 2n   (21 screw)
    if sg == 33:
        if h == 0 and (k + l) % 2 != 0: return False
        if k == 0 and h % 2 != 0: return False
        if h == 0 and k == 0 and l % 2 != 0: return False

    # Pban (#50): b-glide ⊥ a, a-glide ⊥ b, n-glide ⊥ c
    #   0kl: k = 2n, h0l: h = 2n, hk0: h+k = 2n
    if sg == 50:
        if h == 0 and k % 2 != 0: return False
        if k == 0 and h % 2 != 0: return False
        if l == 0 and (h + k) % 2 != 0: return False

    # Pmmn (#59): n-glide ⊥ c
    #   hk0: h+k = 2n
    if sg == 59:
        if l == 0 and (h + k) % 2 != 0: return False

    # Pbcm (#57): b-glide ⊥ a, c-glide ⊥ b, m ⊥ c + 21 screw
    #   0kl: k = 2n, h0l: l = 2n, 00l: l = 2n
    if sg == 57:
        if h == 0 and k % 2 != 0: return False
        if k == 0 and l % 2 != 0: return False
        if h == 0 and k == 0 and l % 2 != 0: return False

    # ── Hexagonal / Trigonal ────────────────────────────────────────────────

    # P63/mmc (#194, e.g. Mo2C, beta-W2C, graphite, WC):
    # (00l): l must be even (63 screw)
    if sg == 194:
        if h == 0 and k == 0 and l % 2 != 0:
            return False

    # P6/mmm (#191), P63/m (#176), P63 (#173), P6322 (#182), P63mc (#186)
    if sg in (173, 176, 182, 186, 191):
        if h == 0 and k == 0 and l % 2 != 0:
            return False

    # P-6m2 (#187, alpha-WC, hexagonal): no screw absences beyond lattice
    # P-43m (#215), F-43m (#216) — already handled by F centering above

    # ── Cubic ──────────────────────────────────────────────────────────────

    # Pm-3n (#223, beta-W A15, Cr3Si type):
    #   hhl: l = 2n (n-glide ⊥ ⟨110⟩, by cubic symmetry)
    #   00l: l = 2n
    if sg == 223:
        ah, ak, al = abs(h), abs(k), abs(l)
        if ah == ak and al % 2 != 0: return False   # hhl: l=2n
        if ah == al and ak % 2 != 0: return False   # hlh → k=2n
        if ak == al and ah % 2 != 0: return False   # hkk → h=2n
        if h == 0 and k == 0 and l % 2 != 0: return False
        if k == 0 and l == 0 and h % 2 != 0: return False
        if h == 0 and l == 0 and k % 2 != 0: return False

    # Ia-3d (#230, garnets): I-centred + d-glide
    if sg == 230:
        if h == 0 and k == 0 and l % 4 != 0: return False  # 00l: l=4n
        if h == 0 and l == 0 and k % 4 != 0: return False
        if k == 0 and l == 0 and h % 4 != 0: return False

    # Fd-3m (#227, diamond, Si, Ge, spinels): F-centred + d-glide
    # F-centering is handled above (h,k,l all odd or all even).
    # d-glide adds: when h,k,l are ALL EVEN, h+k+l must be 4n.
    # Without this, reflections like (200), (222), (442) appear as
    # spurious zero-intensity ticks.
    if sg == 227:
        if h % 2 == 0 and k % 2 == 0 and l % 2 == 0:
            if (h + k + l) % 4 != 0:
                return False

    # Fd-3 (#203): same d-glide rule as 227
    if sg == 203:
        if h % 2 == 0 and k % 2 == 0 and l % 2 == 0:
            if (h + k + l) % 4 != 0:
                return False

    # Im-3m (#229), Pm-3m (#221): handled by lattice centering above

    # ── Orthorhombic glide-plane absences ────────────────────────────────────
    # Critical for transition metal carbides, oxides, and ceramics.

    # Pbcn (#60, W2C, Mo2C):  b-glide ⊥ a, c-glide ⊥ b, n-glide ⊥ c
    #   0kl: k = 2n   |   h0l: l = 2n   |   hk0: h+k = 2n
    if sg == 60:
        if h == 0 and k % 2 != 0:            return False  # 0kl: k=2n
        if k == 0 and l % 2 != 0:            return False  # h0l: l=2n
        if l == 0 and (h + k) % 2 != 0:      return False  # hk0: h+k=2n
        if h == 0 and k == 0 and l % 2 != 0: return False  # 00l: l=2n
        if h == 0 and l == 0 and k % 2 != 0: return False  # 0k0: k=2n
        if k == 0 and l == 0 and h % 2 != 0: return False  # h00: h=2n

    # Pbca (#61, common for phosphates, silicates):
    #   0kl: k = 2n   |   h0l: l = 2n   |   hk0: h = 2n
    if sg == 61:
        if h == 0 and k % 2 != 0:            return False
        if k == 0 and l % 2 != 0:            return False
        if l == 0 and h % 2 != 0:            return False
        if h == 0 and k == 0 and l % 2 != 0: return False
        if h == 0 and l == 0 and k % 2 != 0: return False
        if k == 0 and l == 0 and h % 2 != 0: return False

    # Pnma (#62, Fe3C cementite, many perovskites):
    #   0kl: k+l = 2n   |   hk0: h = 2n   |   h0l: (no condition)
    if sg == 62:
        if h == 0 and (k + l) % 2 != 0:      return False  # 0kl: k+l=2n
        if l == 0 and h % 2 != 0:             return False  # hk0: h=2n
        if h == 0 and k == 0 and l % 2 != 0:  return False  # 00l: l=2n
        if h == 0 and l == 0 and k % 2 != 0:  return False  # 0k0: k=2n
        if k == 0 and l == 0 and h % 2 != 0:  return False  # h00: h=2n

    # P21/c (#14, very common monoclinic):  c-glide ⊥ b, 21 screw along b
    #   h0l: l = 2n   |   0k0: k = 2n
    if sg == 14:
        if k == 0 and l % 2 != 0: return False  # h0l: l=2n
        if h == 0 and l == 0 and k % 2 != 0: return False  # 0k0: k=2n

    # C2/c (#15):  C-centring (h+k even) already handled + c-glide ⊥ b
    #   h0l: l = 2n   |   0k0: k = 2n
    if sg == 15:
        if k == 0 and l % 2 != 0: return False
        if h == 0 and l == 0 and k % 2 != 0: return False

    # P21 21 21 (#19, many organic/MOF crystals):
    #   h00: h = 2n   |   0k0: k = 2n   |   00l: l = 2n
    if sg == 19:
        if k == 0 and l == 0 and h % 2 != 0: return False
        if h == 0 and l == 0 and k % 2 != 0: return False
        if h == 0 and k == 0 and l % 2 != 0: return False

    # R-3m (#166) and R-3c (#167): R-centring (obverse) −h+k+l = 3n
    if sg in (146, 148, 155, 160, 161, 166, 167):
        if (-h + k + l) % 3 != 0: return False
    # R-3c also has: h0l with l=2n, 00l with l=6n (screw)
    if sg in (161, 167):
        if k == 0 and l % 2 != 0: return False
        if h == 0 and k == 0 and l % 6 != 0: return False

    # Pa-3 (#205, pyrite structure):  a-glide
    #   0kl: k = 2n   |   h0l: l = 2n (by cubic equivalence, also h00: h=2n)
    if sg == 205:
        if h == 0 and k % 2 != 0: return False
        if k == 0 and l % 2 != 0: return False
        if l == 0 and h % 2 != 0: return False
        if k == 0 and l == 0 and h % 2 != 0: return False
        if h == 0 and l == 0 and k % 2 != 0: return False
        if h == 0 and k == 0 and l % 2 != 0: return False

    return True

# ─────────────────────────────────────────────────────────────────────────────
# ATOMIC SCATTERING FACTORS  (Cromer-Mann 9-parameter coefficients)
# ─────────────────────────────────────────────────────────────────────────────

# Each entry: (a1,b1, a2,b2, a3,b3, a4,b4, c)
# Source: International Tables for Crystallography, Vol. C, Table 6.1.1.4
# f0(s) = c + Σ_i a_i * exp(-b_i * s²)  where s = sin(θ)/λ in Å⁻¹
CROMER_MANN = {
    'H':  (0.489918,20.6593, 0.262003,7.74039, 0.196767,49.5519, 0.049879,2.20159, 0.001305),
    'C':  (2.31000,20.8439, 1.02000,10.2075, 1.58860,0.568700, 0.865000,51.6512, 0.215600),
    'N':  (12.2126,0.005700, 3.13220,9.89330, 2.01250,28.9975, 1.16630,0.582600,-11.529),
    'O':  (3.04850,13.2771, 2.28680,5.70110, 1.54630,0.323900, 0.867000,32.9089, 0.250800),
    'Si': (6.29150,2.43860, 3.03530,32.3337, 1.98910,0.678500, 1.54100,81.6937, 1.14070),
    'P':  (6.43450,1.90670, 4.17910,27.1570, 1.78000,0.526000, 1.49080,68.1645, 1.11490),
    'S':  (6.90530,1.46790, 5.20340,22.2151, 1.43790,0.253600, 1.58630,56.1720, 0.866900),
    'Ti': (9.75950,7.85080, 7.35580,0.500000, 1.69910,35.6338, 1.90210,116.105, 1.28070),
    'V':  (10.2971,6.86570, 7.35110,0.438500, 2.07030,26.8938, 2.05710,102.478, 1.21990),
    'Cr': (10.6406,6.10380, 7.35370,0.392000, 3.32400,20.2626, 1.49220,98.7399, 1.18320),
    'Mn': (11.2819,5.34090, 7.35730,0.343200, 3.01930,17.8674, 2.24410,83.7543, 1.08960),
    'Fe': (11.7695,4.76110, 7.35730,0.307200, 3.52220,15.3535, 2.30450,76.8805, 1.03690),
    'Co': (12.2841,4.27910, 7.34090,0.278400, 4.00340,13.5359, 2.34880,71.1692, 1.01180),
    'Ni': (12.8376,3.87850, 7.29200,0.256500, 4.44380,12.1763, 2.38000,66.3421, 1.03410),
    'Cu': (13.3380,3.58280, 7.16760,0.247000, 5.61580,11.3966, 1.67350,64.8126, 1.19100),
    'Zn': (14.0743,3.26550, 7.03180,0.233300, 5.16520,10.3163, 2.41000,58.7097, 1.30410),
    'Zr': (17.8765,1.27618, 10.9480,11.9160, 5.41732,0.117622, 3.65721,87.6627, 2.06929),
    'Nb': (17.6142,1.18865, 12.0144,11.7660, 4.04183,0.204785, 3.53346,69.7957, 3.75591),
    'Mo': (3.70250,0.277200, 17.2356,1.09580, 12.8876,11.0040, 3.74290,61.6584, 4.38750),
    'Ru': (19.2674,0.808520, 12.9182,8.43467, 4.86337,24.7997, 1.56756,94.2928, 5.37874),
    'Pd': (19.3319,0.698655, 15.5017,7.98929, 5.29537,25.2052, 0.605844,76.8986, 5.26593),
    'Ag': (19.2808,0.644600, 16.6885,7.47260, 4.80450,24.6605, 1.04630,99.8156, 5.17900),
    'Sn': (19.1889,5.83030, 19.1005,0.503100, 4.45850,26.8909, 2.46630,83.9571, 4.78210),
    'Hf': (29.1440,1.83260, 15.1726,6.70840, 14.7586,0.321800, 4.30013,25.8449,-6.53900),
    'Ta': (29.2024,1.77330, 15.2293,6.43530, 14.5135,0.295100, 4.76492,23.8132,-6.07300),
    'W':  (29.0818,1.72029, 15.4300,6.59432, 14.4327,0.321703, 5.11982,25.2017,-2.98400),
    'Re': (28.7621,1.67190, 15.7189,6.75390, 14.5564,0.334200, 5.44174,24.9604,-2.30400),
    'Os': (28.1894,1.62903, 16.1550,6.83390, 14.9305,0.344000, 5.67589,25.4571,-1.63600),
    'Pt': (27.0059,1.51293, 17.7639,8.81174, 15.7131,0.424593, 5.78370,38.6103, 11.6883),
    'Au': (16.8819,0.461100, 18.5913,8.62160, 25.5582,1.48260, 5.86000,36.3956, 12.0658),
    'Pb': (31.0617,0.690200, 13.0637,2.35760, 18.4420,8.61800, 5.96960,47.2579, 13.4118),
    'Al': (6.42020,3.03870, 1.90020,0.742600, 1.59360,31.5472, 1.96460,85.0886, 1.11510),
    'B':  (2.05450,23.2185, 1.33260,1.02100, 1.09790,60.3498, 0.706800,0.140300,-0.19320),
}


def atomic_scattering_factor(element, s):
    """
    Compute f0 for a given element at s = sin(θ)/λ (in Å⁻¹).
    Falls back to Z/4 approximation if element not in table.
    """
    el = element.strip().capitalize()
    params = CROMER_MANN.get(el)
    if params is None:
        # Rough fallback: use number of electrons / 4
        # (very rough but better than nothing)
        return 10.0
    a1,b1, a2,b2, a3,b3, a4,b4, c = params
    s2 = s * s
    return (c + a1*math.exp(-b1*s2) + a2*math.exp(-b2*s2)
              + a3*math.exp(-b3*s2) + a4*math.exp(-b4*s2))


# ─────────────────────────────────────────────────────────────────────────────
# MOLAR MASS FROM CHEMICAL FORMULA
# ─────────────────────────────────────────────────────────────────────────────

_ATOMIC_MASS = {
    'H':1.008,'He':4.003,'Li':6.941,'Be':9.012,'B':10.81,'C':12.01,'N':14.01,
    'O':16.00,'F':19.00,'Ne':20.18,'Na':22.99,'Mg':24.31,'Al':26.98,'Si':28.09,
    'P':30.97,'S':32.07,'Cl':35.45,'Ar':39.95,'K':39.10,'Ca':40.08,'Sc':44.96,
    'Ti':47.87,'V':50.94,'Cr':52.00,'Mn':54.94,'Fe':55.85,'Co':58.93,'Ni':58.69,
    'Cu':63.55,'Zn':65.38,'Ga':69.72,'Ge':72.63,'As':74.92,'Se':78.97,'Br':79.90,
    'Rb':85.47,'Sr':87.62,'Y':88.91,'Zr':91.22,'Nb':92.91,'Mo':95.95,'Ru':101.1,
    'Rh':102.9,'Pd':106.4,'Ag':107.9,'Cd':112.4,'In':114.8,'Sn':118.7,'Sb':121.8,
    'Te':127.6,'I':126.9,'Cs':132.9,'Ba':137.3,'La':138.9,'Ce':140.1,'Pr':140.9,
    'Nd':144.2,'Sm':150.4,'Eu':152.0,'Gd':157.3,'Tb':158.9,'Dy':162.5,'Ho':164.9,
    'Er':167.3,'Tm':168.9,'Yb':173.0,'Lu':175.0,'Hf':178.5,'Ta':180.9,'W':183.8,
    'Re':186.2,'Os':190.2,'Ir':192.2,'Pt':195.1,'Au':197.0,'Hg':200.6,'Tl':204.4,
    'Pb':207.2,'Bi':209.0,'Th':232.0,'U':238.0,
}


def molar_mass_from_formula(formula):
    """
    Parse a chemical formula string and return molar mass in g/mol.
    Handles CIF-style formulas like 'Mo2 C1', 'Fe3 O4', 'Si O2'.
    Returns None if formula cannot be parsed.
    """
    if not formula:
        return None
    total = 0.0
    # Match element symbol + optional count (integer or decimal)
    for el, count in re.findall(r'([A-Z][a-z]?)\s*(\d*\.?\d*)', formula):
        mass = _ATOMIC_MASS.get(el)
        if mass is None:
            return None
        n = float(count) if count else 1.0
        total += mass * n
    return total if total > 0 else None


def structure_factor_sq(h, k, l, sites, sin_theta_over_lambda):
    """
    Compute |F(hkl)|² from a list of atom sites.

    sites: list of (element, x, y, z, occupancy)
           where x,y,z are fractional coordinates.
    sin_theta_over_lambda: sin(θ)/λ in Å⁻¹

    Returns |F|² (float).
    """
    F_real = 0.0
    F_imag = 0.0
    s = sin_theta_over_lambda
    for (el, x, y, z, occ) in sites:
        f0 = atomic_scattering_factor(el, s)
        phase = 2.0 * math.pi * (h*x + k*y + l*z)
        F_real += occ * f0 * math.cos(phase)
        F_imag += occ * f0 * math.sin(phase)
    return F_real * F_real + F_imag * F_imag


def structure_factor_sq_dw(h, k, l, sites, sin_theta_over_lambda, B_iso_map=None):
    """
    Compute |F(hkl)|² with isotropic Debye-Waller (thermal) factors.

    sites: list of (element, x, y, z, occupancy)
    sin_theta_over_lambda: sin(θ)/λ in Å⁻¹
    B_iso_map: dict {element: B_iso} in Å². If None, B_iso=0 (no thermal damping).
               Typical values: 0.3–1.0 Å² for metals, 0.5–2.0 Å² for oxides.

    The Debye-Waller factor is: T = exp(-B_iso × (sin θ / λ)²)
    This damps high-angle reflections, which is physically correct —
    thermal vibrations blur the electron density and reduce coherent scattering.

    Returns |F|² (float).
    """
    F_real = 0.0
    F_imag = 0.0
    s = sin_theta_over_lambda
    s2 = s * s
    for (el, x, y, z, occ) in sites:
        f0 = atomic_scattering_factor(el, s)
        # Debye-Waller factor
        B = 0.0
        if B_iso_map is not None:
            B = B_iso_map.get(el, B_iso_map.get('_all', 0.0))
        dw = math.exp(-B * s2)
        phase = 2.0 * math.pi * (h*x + k*y + l*z)
        F_real += occ * f0 * dw * math.cos(phase)
        F_imag += occ * f0 * dw * math.sin(phase)
    return F_real * F_real + F_imag * F_imag


# ─────────────────────────────────────────────────────────────────────────────
# PATTERN GENERATION
# ─────────────────────────────────────────────────────────────────────────────

def _expand_sites_by_symmetry(sites, spacegroup_number, a, b, c, al, be, ga):
    """Expand asymmetric-unit sites to full unit cell using pymatgen.

    Returns a list of (element, x, y, z, occupancy) for all atoms in the
    conventional unit cell, or None if expansion fails.
    """
    try:
        from pymatgen.core import Structure, Lattice
        lat = Lattice.from_parameters(a, b, c, al, be, ga)
        species = [s[0] for s in sites]
        coords = [[s[1], s[2], s[3]] for s in sites]
        struct = Structure.from_spacegroup(
            spacegroup_number, lat, species, coords, tol=0.01)
        full_sites = []
        for site in struct:
            fc = site.frac_coords % 1.0
            el = str(site.specie)
            full_sites.append((el, float(fc[0]), float(fc[1]),
                               float(fc[2]), 1.0))
        return full_sites
    except Exception:
        return None


def generate_reflections(a, b, c, al, be, ga, system, spacegroup_number,
                          wavelength, two_theta_min, two_theta_max, hkl_max=10,
                          sites=None, site_policy='auto'):
    """
    Generate list of [two_theta, d, (h,k,l), intensity_weight] for all
    allowed reflections in the given 2θ range.

    site_policy controls how sites are used for structure factor F²:
      'auto' (default) — expand asymmetric unit by symmetry before F²
      'legacy_direct_sites' — use sites directly (old behavior, for
          phases like MP W2C Pbcn where expansion corrupts coordinates)
      'structure_expanded' — always expand (same as auto)

    If `sites` is provided (list of (element, x, y, z, occupancy) tuples),
    the intensity weight is:  multiplicity × |F(hkl)|² × LP
    This correctly produces zero-intensity for structure-factor-extinct peaks.

    If `sites` is None (no atom positions available), the intensity weight
    is just the multiplicity (backward-compatible behaviour).

    Multiplicity is determined by counting every symmetry-equivalent
    (h,k,l) in the full sphere that passes the systematic-absence filter
    and maps to the same d-spacing.
    """
    # ── Expand asymmetric-unit sites to full unit cell ──────────────────
    # structure_factor_sq needs ALL atoms in the unit cell, not just the
    # asymmetric unit.  CIF sites from parse_cif are typically the
    # asymmetric unit, so we expand by space-group symmetry here.
    #
    # site_policy controls this:
    #   'auto'/'structure_expanded' — expand (default, needed for Si etc.)
    #   'legacy_direct_sites' — use sites directly (MP W2C Pbcn compat)
    full_sites = None
    direct_full_cell_sites = (sites is not None
                              and site_policy == 'legacy_direct_sites')
    if sites is not None:
        if direct_full_cell_sites:
            # Use sites as-is. This is used for P1/full-cell MP CIFs where
            # atom positions already encode systematic extinctions in the
            # source setting; applying a separate SG absence table can remove
            # real reflections when that setting differs from the ITA standard.
            full_sites = sites
        elif spacegroup_number > 1:
            full_sites = _expand_sites_by_symmetry(
                sites, spacegroup_number, a, b, c, al, be, ga)
            if full_sites is None:
                full_sites = None  # expansion failed, skip F²

    # Use full_sites for F² if available, else skip F² calculation
    f2_sites = full_sites if full_sites is not None else None

    # Default B_iso for Debye-Waller damping in the theoretical pattern.
    # Without this, high-angle peaks are too intense because thermal
    # motion is ignored.  0.5 Å² is a reasonable room-temperature default
    # for most inorganic materials.
    _DEFAULT_B_ISO = 0.5
    _b_iso_map = {'_all': _DEFAULT_B_ISO}

    seen_d = {}   # d_key → [two_theta, d, (h,k,l), count, weight]

    for h in range(-hkl_max, hkl_max + 1):
        for k in range(-hkl_max, hkl_max + 1):
            for l in range(-hkl_max, hkl_max + 1):
                if h == 0 and k == 0 and l == 0:
                    continue

                if (not direct_full_cell_sites
                        and not is_allowed(h, k, l, spacegroup_number)):
                    continue

                d = d_spacing(h, k, l, a, b, c, al, be, ga, system)
                if d is None or d <= 0:
                    continue

                sin_theta = wavelength / (2 * d)
                if sin_theta > 1 or sin_theta < 0:
                    continue

                two_theta = 2 * math.degrees(math.asin(sin_theta))
                if two_theta < two_theta_min or two_theta > two_theta_max:
                    continue

                d_key = round(d, 4)
                if d_key in seen_d:
                    seen_d[d_key][3] += 1        # count this equivalent
                else:
                    # Compute |F|² with Debye-Waller damping, plus LP factor,
                    # so the theoretical stick intensities match a measured
                    # powder pattern visually.
                    F2_LP = None
                    if f2_sites is not None:
                        s = sin_theta / wavelength
                        F2 = structure_factor_sq_dw(h, k, l, f2_sites, s,
                                                     _b_iso_map)
                        # Lorentz-polarisation factor
                        theta_r = math.radians(two_theta / 2)
                        cos2t = math.cos(math.radians(two_theta))
                        sin_t = math.sin(theta_r)
                        cos_t = math.cos(theta_r)
                        LP = ((1 + cos2t**2) / (sin_t**2 * cos_t)
                              if sin_t > 0 and cos_t > 0 else 1.0)
                        F2_LP = F2 * LP
                    seen_d[d_key] = [two_theta, d, (abs(h), abs(k), abs(l)),
                                     1, F2_LP]

    # Build final list with intensity weights.
    entries = sorted(seen_d.values(), key=lambda x: x[0])

    max_w = 0.0
    if f2_sites is not None:
        for entry in entries:
            w = entry[4]
            if w is not None and w > max_w:
                max_w = w
    rel_threshold = max_w * 1e-2  # 1% of strongest reflection

    reflections = []
    for entry in entries:
        two_theta, d, hkl, mult, F2_LP = entry
        if f2_sites is not None and F2_LP is not None:
            # Skip reflections with negligible intensity
            if F2_LP < max(1e-4, rel_threshold):
                continue
            weight = mult * F2_LP
        else:
            # No F² available — use multiplicity only
            weight = mult
        reflections.append([two_theta, d, hkl, weight])

    return reflections  # [(two_theta, d, (h,k,l), intensity_weight), ...]


def generate_reflections_rietveld(a, b, c, al, be, ga, system, spacegroup_number,
                                   wavelength, two_theta_min, two_theta_max,
                                   sites, hkl_max=12, site_policy='auto'):
    """
    Generate reflection list for Rietveld refinement.

    Unlike generate_reflections(), this returns per-reflection components
    that allow recomputing intensities when B_iso changes without
    re-enumerating all hkl:

    site_policy: 'auto', 'legacy_direct_sites', or 'structure_expanded'
        (see generate_reflections docstring).

    Returns list of dicts:
        {'two_theta', 'd', 'hkl', 'mult', 'sin_theta_over_lambda',
         'h_rep', 'k_rep', 'l_rep'}

    The representative (h,k,l) is stored so structure_factor_sq_dw() can
    be called with varying B_iso during refinement.
    """
    # Expand asymmetric-unit sites to full unit cell for correct F²
    full_sites = sites
    if sites and spacegroup_number > 1:
        if site_policy == 'legacy_direct_sites':
            full_sites = sites
        else:
            expanded = _expand_sites_by_symmetry(
                sites, spacegroup_number, a, b, c, al, be, ga)
            if expanded is not None:
                full_sites = expanded

    seen_d = {}

    for h in range(-hkl_max, hkl_max + 1):
        for k in range(-hkl_max, hkl_max + 1):
            for l in range(-hkl_max, hkl_max + 1):
                if h == 0 and k == 0 and l == 0:
                    continue
                if not is_allowed(h, k, l, spacegroup_number):
                    continue

                d = d_spacing(h, k, l, a, b, c, al, be, ga, system)
                if d is None or d <= 0:
                    continue

                sin_theta = wavelength / (2 * d)
                if sin_theta > 1 or sin_theta < 0:
                    continue

                two_theta = 2 * math.degrees(math.asin(sin_theta))
                if two_theta < two_theta_min or two_theta > two_theta_max:
                    continue

                d_key = round(d, 4)
                if d_key in seen_d:
                    seen_d[d_key]['mult'] += 1
                else:
                    # Check if this reflection has nonzero F²
                    s = sin_theta / wavelength
                    F2_test = structure_factor_sq(h, k, l, full_sites, s)
                    if F2_test < 1e-4:
                        seen_d[d_key] = {'skip': True, 'mult': 1}
                        continue
                    seen_d[d_key] = {
                        'two_theta': two_theta,
                        'd': d,
                        'hkl': (abs(h), abs(k), abs(l)),
                        'mult': 1,
                        'sin_theta_over_lambda': s,
                        'h_rep': h, 'k_rep': k, 'l_rep': l,
                        'skip': False,
                    }

    refs = [v for v in sorted(seen_d.values(), key=lambda x: x.get('two_theta', 999))
            if not v.get('skip', False)]
    return refs


def compute_rietveld_intensities(refs, sites, B_iso_map=None,
                                  spacegroup_number=None,
                                  a=None, b=None, c=None,
                                  al=None, be=None, ga=None):
    """
    Compute per-reflection intensities for a Rietveld model:
      I_hkl = mult × |F(hkl)|² × LP

    refs: output of generate_reflections_rietveld()
    sites: atom site list (asymmetric unit or full cell)
    B_iso_map: {element: B_iso} or {'_all': B_iso}
    spacegroup_number, a, b, c, al, be, ga: if provided, sites are
        expanded from asymmetric unit to full cell before computing F².
        This ensures consistency with generate_reflections_rietveld().

    Returns np.array of intensities, same length as refs.
    """
    # Expand asymmetric-unit sites to full cell if SG info provided
    full_sites = sites
    if (sites and spacegroup_number is not None and spacegroup_number > 1
            and a is not None):
        expanded = _expand_sites_by_symmetry(
            sites, spacegroup_number,
            a, b or a, c or a,
            al or 90.0, be or 90.0, ga or 90.0)
        if expanded is not None:
            full_sites = expanded

    intensities = np.zeros(len(refs))
    for i, ref in enumerate(refs):
        s = ref['sin_theta_over_lambda']
        h, k, l = ref['h_rep'], ref['k_rep'], ref['l_rep']
        F2 = structure_factor_sq_dw(h, k, l, full_sites, s, B_iso_map)
        # Lorentz-polarisation factor
        tt = ref['two_theta']
        theta_r = math.radians(tt / 2)
        cos2t = math.cos(math.radians(tt))
        sin_t = math.sin(theta_r)
        cos_t = math.cos(theta_r)
        LP = (1 + cos2t**2) / (sin_t**2 * cos_t) if sin_t > 0 and cos_t > 0 else 1.0
        intensities[i] = ref['mult'] * F2 * LP
    return intensities


# ─────────────────────────────────────────────────────────────────────────────
# PEAK SHAPE: THOMPSON-COX-HASTINGS PSEUDO-VOIGT
# ─────────────────────────────────────────────────────────────────────────────

def caglioti_fwhm(two_theta_deg, U, V, W):
    """
    Caglioti equation: FWHM² = U·tan²θ + V·tanθ + W
    Returns FWHM in degrees.
    """
    theta = math.radians(two_theta_deg / 2)
    tan_t = math.tan(theta)
    fwhm2 = U * tan_t**2 + V * tan_t + W
    if fwhm2 <= 0:
        return 0.01
    return math.sqrt(fwhm2)


def tch_fwhm_eta(two_theta_deg, U, V, W, X, Y):
    """
    Thompson-Cox-Hastings pseudo-Voigt profile parameters.

    Gaussian FWHM:  H_G² = U·tan²θ + V·tanθ + W
    Lorentzian FWHM: H_L = X·tanθ + Y/cosθ

    Where:
      U = Gaussian micro-strain broadening
      V, W = instrumental Gaussian parameters
      X = Lorentzian micro-strain broadening
      Y = Lorentzian size broadening (Scherrer term: Y = Kλ/L in degrees)

    Returns (total_FWHM, eta) using the TCH mixing approximation.
    """
    theta = math.radians(two_theta_deg / 2)
    tan_t = math.tan(theta)
    cos_t = math.cos(theta)

    # Gaussian component
    hg2 = U * tan_t**2 + V * tan_t + W
    H_G = math.sqrt(max(hg2, 1e-8))

    # Lorentzian component
    H_L = max(X * tan_t + Y / max(cos_t, 1e-8), 1e-6)

    # TCH approximation for total FWHM  (Thompson, Cox & Hastings 1987)
    hg5 = H_G**5
    hl5 = H_L**5
    H5 = (hg5
           + 2.69269 * H_G**4 * H_L
           + 2.42843 * H_G**3 * H_L**2
           + 4.47163 * H_G**2 * H_L**3
           + 0.07842 * H_G    * H_L**4
           + hl5)
    H = H5 ** 0.2  # fifth root

    # Mixing parameter η  (fraction Lorentzian)
    q = H_L / max(H, 1e-8)
    eta = max(0.0, min(1.0,
              1.36603 * q - 0.47719 * q**2 + 0.11116 * q**3))

    return max(H, 0.005), eta


def pseudo_voigt(x, x0, fwhm, eta):
    """
    Pseudo-Voigt profile: η·Lorentzian + (1-η)·Gaussian
    eta = 0 → pure Gaussian, eta = 1 → pure Lorentzian
    """
    eta = max(0.0, min(1.0, eta))
    sigma_g = fwhm / (2 * math.sqrt(2 * math.log(2)))
    gamma_l = fwhm / 2

    dx = x - x0
    gauss = np.exp(-0.5 * (dx / sigma_g)**2)
    lor   = 1.0 / (1.0 + (dx / gamma_l)**2)
    return eta * lor + (1 - eta) * gauss

def compute_phase_pattern(two_theta_array, reflections, scale,
                           U, V, W, eta=None, two_theta_zero=0.0,
                           X=0.0, Y=0.0):
    """
    Compute the simulated intensity pattern for one phase.

    Profile model:
      If X or Y are non-zero, uses TCH pseudo-Voigt where eta is computed
      per-peak from the Gaussian (U,V,W) and Lorentzian (X,Y) widths.
      Otherwise falls back to fixed-eta Caglioti pseudo-Voigt.

    reflections: output of generate_reflections()
    Returns array of intensities same shape as two_theta_array.
    """
    use_tch = (X != 0.0 or Y != 0.0)
    pattern = np.zeros_like(two_theta_array)
    tt_shifted = two_theta_array - two_theta_zero

    for tt_peak, d, hkl, mult in reflections:
        if use_tch:
            fwhm, eta_pk = tch_fwhm_eta(tt_peak, U, V, W, X, Y)
        else:
            fwhm = max(caglioti_fwhm(tt_peak, U, V, W), 0.005)
            eta_pk = eta if eta is not None else 0.5
        # Only compute within ±20*fwhm of peak centre (efficiency)
        window = 20 * fwhm
        mask = np.abs(tt_shifted - tt_peak) < window
        if not mask.any():
            continue
        profile = pseudo_voigt(tt_shifted[mask], tt_peak, fwhm, eta_pk)
        # Lorentz-polarisation factor
        theta_r = math.radians(tt_peak / 2)
        cos2t    = math.cos(math.radians(tt_peak))
        Lp = (1 + cos2t**2) / (math.sin(theta_r)**2 * math.cos(theta_r))
        pattern[mask] += scale * mult * Lp * profile

    return pattern


# ─────────────────────────────────────────────────────────────────────────────
# CHEBYSHEV POLYNOMIAL BACKGROUND
# ─────────────────────────────────────────────────────────────────────────────

def chebyshev_background(two_theta, coeffs, tt_min, tt_max):
    """
    Evaluate Chebyshev polynomial background.
    coeffs: list of polynomial coefficients [c0, c1, c2, ...]
    Normalises 2θ to [-1, 1] for numerical stability.
    """
    x = 2 * (two_theta - tt_min) / (tt_max - tt_min) - 1
    bg = np.zeros_like(two_theta)
    T_prev2 = np.ones_like(x)      # T0
    T_prev1 = x.copy()             # T1
    if len(coeffs) > 0:
        bg += coeffs[0] * T_prev2
    if len(coeffs) > 1:
        bg += coeffs[1] * T_prev1
    for i, c in enumerate(coeffs[2:], 2):
        T_curr = 2 * x * T_prev1 - T_prev2
        bg += c * T_curr
        T_prev2, T_prev1 = T_prev1, T_curr
    return bg


# ─────────────────────────────────────────────────────────────────────────────
# SCHERRER CRYSTALLITE SIZE
# ─────────────────────────────────────────────────────────────────────────────

def scherrer_size(fwhm_deg, two_theta_deg, wavelength, K=0.94):
    """
    Scherrer equation: D = K·λ / (β·cos θ)
    fwhm_deg: FWHM in degrees (instrumental broadening NOT subtracted here)
    Returns crystallite size D in Å.
    """
    beta_rad = math.radians(fwhm_deg)
    theta_rad = math.radians(two_theta_deg / 2)
    if beta_rad <= 0 or math.cos(theta_rad) <= 0:
        return None
    return K * wavelength / (beta_rad * math.cos(theta_rad))


def size_from_Y(Y_deg, wavelength, K=0.94):
    """
    Extract crystallite size directly from the TCH Lorentzian size parameter Y.

    In the TCH model, Y/cosθ is the Lorentzian size contribution to FWHM,
    which is exactly the Scherrer term:  Y = Kλ / (L · π/180)
    Solving: L = Kλ / (Y · π/180)

    Returns crystallite size L in Å, or None if Y ≤ 0.
    """
    if Y_deg <= 0:
        return None
    return K * wavelength / math.radians(Y_deg)


# ─────────────────────────────────────────────────────────────────────────────
# SYMMETRY EXPANSION
# ─────────────────────────────────────────────────────────────────────────────

def expand_sites_from_cif(cif_text):
    """Expand CIF asymmetric-unit sites to the full unit cell using pymatgen.

    CIF files from COD and many databases list only the *asymmetric unit*.
    Computing |F(hkl)|² from those sites alone gives incorrect structure
    factors (e.g. F(110)≠0 for A15-W when it should be 0).

    Returns a list of (element, x, y, z, occupancy) tuples covering the
    full conventional unit cell, or *None* if expansion fails.
    """
    if not cif_text:
        return None
    try:
        from pymatgen.io.cif import CifParser
        try:
            parser = CifParser.from_str(cif_text, occupancy_tolerance=100.0)
        except (AttributeError, TypeError):
            import tempfile as _tf, os as _os
            _fd, _tmp = _tf.mkstemp(suffix='.cif')
            with _os.fdopen(_fd, 'w') as _f:
                _f.write(cif_text)
            parser = CifParser(_tmp, occupancy_tolerance=100.0)
            _os.unlink(_tmp)
        structs = parser.parse_structures(primitive=False)
        if structs:
            sites = []
            for site in structs[0]:
                frac = site.frac_coords % 1.0
                el = str(site.specie)
                occ = 1.0
                if hasattr(site, 'properties') and 'occupancy' in site.properties:
                    occ = float(site.properties['occupancy'])
                sites.append((el, float(frac[0]), float(frac[1]),
                              float(frac[2]), occ))
            if sites:
                return sites
    except BaseException:
        pass

    # ── Fallback: expand using built-in symmetry operations ─────────────
    # When pymatgen is unavailable, apply space group symmetry operations
    # to the asymmetric unit parsed from the CIF.
    try:
        parsed = parse_cif(cif_text)
        asym_sites = parsed.get('sites') or []
        sg = parsed.get('spacegroup_number', 1)
        if asym_sites and sg > 1:
            expanded = _expand_by_symmetry(asym_sites, sg)
            if expanded and len(expanded) > len(asym_sites):
                return expanded
    except Exception:
        pass

    return None


# ─────────────────────────────────────────────────────────────────────────────
# BUILT-IN SYMMETRY EXPANSION (fallback when pymatgen is unavailable)
# ─────────────────────────────────────────────────────────────────────────────

# Symmetry operations as (rotation_matrix_rows, translation_vector).
# Each operation: ((r00,r01,r02), (r10,r11,r12), (r20,r21,r22)), (tx,ty,tz)
# Only the most common space groups are included; others return None.

def _apply_symop(site, op):
    """Apply a symmetry operation to a fractional coordinate site."""
    el, x, y, z, occ = site
    rot, trans = op
    xn = rot[0][0]*x + rot[0][1]*y + rot[0][2]*z + trans[0]
    yn = rot[1][0]*x + rot[1][1]*y + rot[1][2]*z + trans[1]
    zn = rot[2][0]*x + rot[2][1]*y + rot[2][2]*z + trans[2]
    return (el, xn % 1.0, yn % 1.0, zn % 1.0, occ)


def _expand_by_symmetry(asym_sites, spacegroup_number):
    """Expand asymmetric unit to full unit cell using space group operations.

    Returns list of (element, x, y, z, occupancy) for the full unit cell,
    or None if the space group is not in our built-in table.
    """
    ops = _SG_SYMOPS.get(spacegroup_number)
    if ops is None:
        return None

    all_sites = []
    for site in asym_sites:
        seen = set()
        for op in ops:
            new_site = _apply_symop(site, op)
            # Round to avoid floating-point duplicates
            key = (new_site[0],
                   round(new_site[1], 4),
                   round(new_site[2], 4),
                   round(new_site[3], 4))
            if key not in seen:
                seen.add(key)
                all_sites.append(new_site)
    return all_sites if all_sites else None


# Identity operation (used by all space groups)
_E = ((1,0,0),(0,1,0),(0,0,1)), (0,0,0)

# Symmetry operations for common space groups.
# These are the FULL set of general-position operations from ITC Vol A.

_SG_SYMOPS = {
    # P1 (#1): just identity
    1: [_E],

    # P-1 (#2): identity + inversion
    2: [
        _E,
        (((-1,0,0),(0,-1,0),(0,0,-1)), (0,0,0)),
    ],

    # P21/c (#14): 4 operations
    14: [
        _E,
        (((-1,0,0),(0,1,0),(0,0,-1)), (0,0.5,0.5)),   # 21 screw + c-glide
        (((-1,0,0),(0,-1,0),(0,0,-1)), (0,0,0)),        # inversion
        (((1,0,0),(0,-1,0),(0,0,1)), (0,0.5,0.5)),      # glide
    ],

    # Pbcn (#60): 8 operations  (ITC Vol A, origin choice 1)
    60: [
        _E,
        (((-1,0,0),(0,-1,0),(0,0,1)),  (0.5,0.5,0.5)),  # 2₁(001)
        (((1,0,0),(0,-1,0),(0,0,-1)),  (0.5,0.5,0.5)),   # 2₁(010)
        (((-1,0,0),(0,1,0),(0,0,-1)),  (0,0,0)),          # 2(100)
        (((-1,0,0),(0,-1,0),(0,0,-1)), (0,0,0)),          # inversion
        (((1,0,0),(0,1,0),(0,0,-1)),   (0.5,0.5,0.5)),   # b-glide
        (((-1,0,0),(0,1,0),(0,0,1)),   (0.5,0.5,0.5)),   # c-glide
        (((1,0,0),(0,-1,0),(0,0,1)),   (0,0,0)),          # n-glide
    ],

    # Pbca (#61): 8 operations
    61: [
        _E,
        (((-1,0,0),(0,-1,0),(0,0,1)),  (0.5,0,0.5)),
        (((1,0,0),(0,-1,0),(0,0,-1)),  (0,0.5,0.5)),
        (((-1,0,0),(0,1,0),(0,0,-1)),  (0.5,0.5,0)),
        (((-1,0,0),(0,-1,0),(0,0,-1)), (0,0,0)),
        (((1,0,0),(0,1,0),(0,0,-1)),   (0.5,0,0.5)),
        (((-1,0,0),(0,1,0),(0,0,1)),   (0,0.5,0.5)),
        (((1,0,0),(0,-1,0),(0,0,1)),   (0.5,0.5,0)),
    ],

    # Pnma (#62): 8 operations
    62: [
        _E,
        (((-1,0,0),(0,-1,0),(0,0,1)),  (0.5,0,0.5)),
        (((-1,0,0),(0,1,0),(0,0,-1)),  (0,0.5,0)),
        (((1,0,0),(0,-1,0),(0,0,-1)),  (0.5,0.5,0.5)),
        (((-1,0,0),(0,-1,0),(0,0,-1)), (0,0,0)),
        (((1,0,0),(0,1,0),(0,0,-1)),   (0.5,0,0.5)),
        (((1,0,0),(0,-1,0),(0,0,1)),   (0,0.5,0)),
        (((-1,0,0),(0,1,0),(0,0,1)),   (0.5,0.5,0.5)),
    ],

    # Cmcm (#63): C-centred, 16 positions in full cell
    63: [
        _E,
        (((-1,0,0),(0,-1,0),(0,0,1)),  (0,0,0.5)),
        (((1,0,0),(0,-1,0),(0,0,-1)),  (0,0,0)),
        (((-1,0,0),(0,1,0),(0,0,-1)),  (0,0,0.5)),
        (((-1,0,0),(0,-1,0),(0,0,-1)), (0,0,0)),
        (((1,0,0),(0,1,0),(0,0,-1)),   (0,0,0.5)),
        (((-1,0,0),(0,1,0),(0,0,1)),   (0,0,0)),
        (((1,0,0),(0,-1,0),(0,0,1)),   (0,0,0.5)),
        # C-centering translations of the above 8
        (((1,0,0),(0,1,0),(0,0,1)),    (0.5,0.5,0)),
        (((-1,0,0),(0,-1,0),(0,0,1)),  (0.5,0.5,0.5)),
        (((1,0,0),(0,-1,0),(0,0,-1)),  (0.5,0.5,0)),
        (((-1,0,0),(0,1,0),(0,0,-1)),  (0.5,0.5,0.5)),
        (((-1,0,0),(0,-1,0),(0,0,-1)), (0.5,0.5,0)),
        (((1,0,0),(0,1,0),(0,0,-1)),   (0.5,0.5,0.5)),
        (((-1,0,0),(0,1,0),(0,0,1)),   (0.5,0.5,0)),
        (((1,0,0),(0,-1,0),(0,0,1)),   (0.5,0.5,0.5)),
    ],

    # P63/mmc (#194): 24 operations
    194: [
        _E,
        (((-1,1,0),(-1,0,0),(0,0,1)),  (0,0,0.5)),
        (((0,1,0),(-1,1,0),(0,0,1)),   (0,0,0)),
        (((-1,0,0),(0,-1,0),(0,0,1)),   (0,0,0.5)),
        (((1,-1,0),(1,0,0),(0,0,1)),    (0,0,0)),
        (((0,-1,0),(1,-1,0),(0,0,1)),   (0,0,0.5)),
        (((0,-1,0),(-1,0,0),(0,0,-1)),  (0,0,0.5)),
        (((-1,1,0),(0,1,0),(0,0,-1)),   (0,0,0)),
        (((1,0,0),(1,-1,0),(0,0,-1)),   (0,0,0.5)),
        (((0,1,0),(1,0,0),(0,0,-1)),    (0,0,0)),
        (((1,-1,0),(0,-1,0),(0,0,-1)),  (0,0,0.5)),
        (((-1,0,0),(-1,1,0),(0,0,-1)),  (0,0,0)),
        (((-1,0,0),(0,-1,0),(0,0,-1)),  (0,0,0)),
        (((1,-1,0),(1,0,0),(0,0,-1)),   (0,0,0.5)),
        (((0,-1,0),(1,-1,0),(0,0,-1)),  (0,0,0)),
        (((1,0,0),(0,1,0),(0,0,-1)),    (0,0,0.5)),
        (((-1,1,0),(-1,0,0),(0,0,-1)),  (0,0,0)),
        (((0,1,0),(-1,1,0),(0,0,-1)),   (0,0,0.5)),
        (((0,1,0),(1,0,0),(0,0,1)),     (0,0,0.5)),
        (((1,-1,0),(0,-1,0),(0,0,1)),   (0,0,0)),
        (((-1,0,0),(-1,1,0),(0,0,1)),   (0,0,0.5)),
        (((0,-1,0),(-1,0,0),(0,0,1)),   (0,0,0)),
        (((-1,1,0),(0,1,0),(0,0,1)),    (0,0,0.5)),
        (((1,0,0),(1,-1,0),(0,0,1)),    (0,0,0)),
    ],

    # Im-3m (#229): BCC metals (W, Mo, Fe) — CIFs typically provide
    # all sites already. Rely on is_allowed() I-centering filter.

    # Pm-3n (#223): A15 phases — CIFs typically provide all sites.

    # Fm-3m (#225): FCC metals — CIFs provide full conventional cell.
}


# ─────────────────────────────────────────────────────────────────────────────
# CIF PARSER
# ─────────────────────────────────────────────────────────────────────────────

def parse_cif(cif_text):
    """
    Minimal CIF parser — extracts unit cell, space group, formula,
    and atom site coordinates (fractional).

    Returns dict with keys: a, b, c, alpha, beta, gamma, spacegroup_number,
    spacegroup_name, system, formula, cod_id, sites

    sites: list of (element, x, y, z, occupancy) tuples parsed from the
    _atom_site loop. Empty list if no atom sites found.
    """
    result = {
        'a': None, 'b': None, 'c': None,
        'alpha': 90.0, 'beta': 90.0, 'gamma': 90.0,
        'spacegroup_number': 1,
        'spacegroup_name': 'P 1',
        'system': 'triclinic',
        'formula': '',
        'cod_id': '',
        'sites': [],
        'Z': None,
    }

    def parse_val(s):
        """Strip uncertainty e.g. '3.002(5)' → 3.002"""
        s = s.strip().strip("'\"")
        m = re.match(r'^([0-9\.\-\+eE]+)', s)
        return float(m.group(1)) if m else None

    lines = cif_text.splitlines()

    # ── Pass 1: scalar CIF tags ──────────────────────────────────────────
    for line in lines:
        line = line.strip()
        if line.startswith('_cell_length_a'):
            v = parse_val(line.split()[-1])
            if v: result['a'] = v
        elif line.startswith('_cell_length_b'):
            v = parse_val(line.split()[-1])
            if v: result['b'] = v
        elif line.startswith('_cell_length_c'):
            v = parse_val(line.split()[-1])
            if v: result['c'] = v
        elif line.startswith('_cell_angle_alpha'):
            v = parse_val(line.split()[-1])
            if v: result['alpha'] = v
        elif line.startswith('_cell_angle_beta'):
            v = parse_val(line.split()[-1])
            if v: result['beta'] = v
        elif line.startswith('_cell_angle_gamma'):
            v = parse_val(line.split()[-1])
            if v: result['gamma'] = v
        elif (line.startswith('_symmetry_Int_Tables_number') or
              line.startswith('_space_group_IT_number') or
              line.startswith('_space_group.IT_number')):
            v = parse_val(line.split()[-1])
            if v: result['spacegroup_number'] = int(v)
        elif (line.startswith('_symmetry_space_group_name_H-M') or
              line.startswith('_space_group_name_H-M_alt') or
              line.startswith('_space_group.name_H-M') or
              line.startswith('_space_group_name_H-M ')):
            parts = line.split(None, 1)
            if len(parts) > 1:
                result['spacegroup_name'] = parts[1].strip().strip("'\"")
        elif line.startswith('_chemical_formula_sum'):
            parts = line.split(None, 1)
            if len(parts) > 1:
                result['formula'] = parts[1].strip().strip("'\"")
        elif line.startswith('_cell_formula_units_Z'):
            v = parse_val(line.split()[-1])
            if v: result['Z'] = int(v)

    # ── Pass 2: parse atom_site loop ─────────────────────────────────────
    sites = _parse_atom_site_loop(lines, parse_val)
    if sites:
        result['sites'] = sites

    # Infer crystal system from space group number
    sg = result['spacegroup_number']
    if   1  <= sg <=   2: result['system'] = 'triclinic'
    elif 3  <= sg <=  15: result['system'] = 'monoclinic'
    elif 16 <= sg <=  74: result['system'] = 'orthorhombic'
    elif 75 <= sg <= 142: result['system'] = 'tetragonal'
    elif 143 <= sg <= 167: result['system'] = 'trigonal'
    elif 168 <= sg <= 194: result['system'] = 'hexagonal'
    elif 195 <= sg <= 230: result['system'] = 'cubic'

    return result


def _parse_atom_site_loop(lines, parse_val):
    """
    Extract atom sites from a CIF _atom_site loop block.
    Returns list of (element, x, y, z, occupancy) tuples.

    Handles both _atom_site_type_symbol and _atom_site_label as the
    element identifier, and tolerates varying column orders.
    """
    sites = []
    i = 0
    n = len(lines)

    while i < n:
        line = lines[i].strip()
        # Find a loop_ that contains _atom_site_ tags
        if line.lower() == 'loop_':
            i += 1
            # Collect column headers
            headers = []
            while i < n:
                hline = lines[i].strip()
                if hline.startswith('_atom_site_'):
                    headers.append(hline.lower())
                    i += 1
                else:
                    break
            if not headers:
                continue

            # Identify column indices for fields we need
            col_type  = None  # _atom_site_type_symbol (preferred)
            col_label = None  # _atom_site_label (fallback for element)
            col_x = col_y = col_z = col_occ = None
            for ci, h in enumerate(headers):
                if h == '_atom_site_type_symbol':     col_type  = ci
                elif h == '_atom_site_label':          col_label = ci
                elif h == '_atom_site_fract_x':        col_x     = ci
                elif h == '_atom_site_fract_y':        col_y     = ci
                elif h == '_atom_site_fract_z':        col_z     = ci
                elif h == '_atom_site_occupancy':      col_occ   = ci

            # Must have at least fractional coordinates and some element id
            el_col = col_type if col_type is not None else col_label
            if el_col is None or col_x is None or col_y is None or col_z is None:
                continue

            # Read data rows
            while i < n:
                dline = lines[i].strip()
                if not dline or dline.startswith('_') or dline.lower() == 'loop_' or dline.startswith('#'):
                    break
                parts = dline.split()
                if len(parts) <= max(el_col, col_x, col_y, col_z):
                    break
                try:
                    el_raw = parts[el_col].strip("'\"")
                    # Extract element symbol: strip trailing digits/charges
                    el = re.match(r'([A-Z][a-z]?)', el_raw)
                    if not el:
                        i += 1
                        continue
                    element = el.group(1)
                    x = parse_val(parts[col_x])
                    y = parse_val(parts[col_y])
                    z = parse_val(parts[col_z])
                    occ = parse_val(parts[col_occ]) if col_occ is not None and col_occ < len(parts) else 1.0
                    if x is not None and y is not None and z is not None:
                        sites.append((element, x, y, z, occ if occ is not None else 1.0))
                except (IndexError, ValueError):
                    pass
                i += 1
            continue
        i += 1

    return sites


# ─────────────────────────────────────────────────────────────────────────────
# GOODNESS-OF-FIT STATISTICS
# ─────────────────────────────────────────────────────────────────────────────

def compute_fit_statistics(y_obs, y_calc, y_weights, n_params):
    """
    Compute Rwp, Rp, chi-squared, and GoF.
    y_weights = 1/sigma^2
    """
    diff     = y_obs - y_calc
    Rwp      = math.sqrt(np.sum(y_weights * diff**2) / np.sum(y_weights * y_obs**2))
    Rp       = np.sum(np.abs(diff)) / np.sum(np.abs(y_obs))
    chi2     = np.sum(y_weights * diff**2) / max(len(y_obs) - n_params, 1)
    GoF      = math.sqrt(chi2)
    return {'Rwp': round(Rwp*100,2), 'Rp': round(Rp*100,2),
            'chi2': round(chi2,3), 'GoF': round(GoF,3)}

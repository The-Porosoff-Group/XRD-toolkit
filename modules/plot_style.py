"""Shared Arial typography and display-only scientific labels for PNG exports."""

from functools import lru_cache, wraps
import math
import re
from threading import RLock
import warnings

from PIL import Image, ImageDraw, ImageFont


@lru_cache(maxsize=1)
def plot_font_family():
    """Prefer installed Arial; never redistribute proprietary font files."""
    from matplotlib import font_manager

    for family in ('Arial', 'Liberation Sans', 'DejaVu Sans'):
        try:
            font_manager.findfont(family, fallback_to_default=False)
        except ValueError:
            continue
        if family != 'Arial':
            warnings.warn(f'Arial is not installed; figure exports use {family}.',
                          RuntimeWarning, stacklevel=2)
        return family
    raise RuntimeError('No supported sans-serif plot font is installed.')


@lru_cache(maxsize=128)
def load_plot_font(size, bold=False):
    from matplotlib import font_manager

    properties = font_manager.FontProperties(
        family=plot_font_family(), weight='bold' if bold else 'normal')
    path = font_manager.findfont(properties, fallback_to_default=False)
    return ImageFont.truetype(path, max(7, round(size)))


_MATPLOTLIB_LOCK = RLock()


def arial_plot(function):
    """Scope font settings to an export and serialize Matplotlib rendering."""
    @wraps(function)
    def render(*args, **kwargs):
        import matplotlib

        family = plot_font_family()
        settings = {
            'font.family': family, 'mathtext.fontset': 'custom',
            'mathtext.rm': family, 'mathtext.it': f'{family}:italic',
            'mathtext.bf': f'{family}:bold',
            'mathtext.bfit': f'{family}:italic:bold',
            'mathtext.sf': family, 'mathtext.tt': family,
            'mathtext.cal': family, 'mathtext.fallback': 'stixsans',
            'text.usetex': False,
        }
        with _MATPLOTLIB_LOCK, matplotlib.rc_context(settings):
            return function(*args, **kwargs)
    return render


_SUB = '₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎ₓ'
_SUP = '⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾'
_TO_SUB = str.maketrans('0123456789-x', _SUB[:10] + '₋ₓ')
_FROM_SUB = str.maketrans(_SUB, '0123456789+-=()x')
_FROM_SUP = str.maketrans(_SUP, '0123456789+-=()')
_ELEMENTS = set(('H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V '
                 'Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc '
                 'Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu '
                 'Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi '
                 'Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr '
                 'Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og').split())
_FORMULA = re.compile(
    r'(?<![A-Za-z0-9])(?P<prefix>(?:n|i|c[12]|t[12])-?)?'
    r'(?P<formula>(?:[A-Z][a-z]?(?:\d+(?:\.\d+)?(?:-x)?)?)+)(?![A-Za-z0-9])')
_SCRIPT = re.compile(
    rf'[{_SUB}]+(?:\.[{_SUB}]+)*|[{_SUP}]+|[0-9][\u0304\u0305]')
_EXPLICIT_SCRIPT = re.compile(r'([_^])\{([^{}\r\n]+)\}')
_LABEL_WORD = re.compile(r'(?:[_^]\{[^{}\r\n]+\}|\S)+')


def scientific_unicode(text):
    """Recognize formula counts, preserving run numbers, temperatures and dates.

    Single-element labels are deliberately limited to common molecular gases
    and carbon-number groups. A sample identifier such as W123 stays literal.
    Existing Unicode sub/superscripts also remain valid explicit notation.
    """
    def formula(match):
        value = match['formula']
        elements = re.findall(r'[A-Z][a-z]?', value)
        if (not set(elements) <= _ELEMENTS
                or not re.search(r'\d', value)
                or (len(elements) == 1 and not re.fullmatch(
                    r'(?:[CHNOS]\d+|Cl2|Br2|I2|F2)', value))):
            return match.group(0)
        return (match['prefix'] or '') + value.translate(_TO_SUB)

    return _FORMULA.sub(formula, str(text or ''))


def scientific_runs(text, *, recognize_formulas=True):
    """Typeset formulas, Unicode scripts, and explicit _{text}/^{text} codes.

    Braces are required so ordinary underscores in sample IDs stay literal.
    Explicit contents are plain text, including letters and decimal counts;
    they do not require a TeX installation or change the underlying metadata.
    """
    runs = []

    def append(value, style):
        if not value:
            return
        if runs and runs[-1][1] == style and style != 'bar':
            runs[-1] = (runs[-1][0] + value, style)
        else:
            runs.append((value, style))

    def plain(value):
        value = scientific_unicode(value) if recognize_formulas else value
        end = 0
        for match in _SCRIPT.finditer(value):
            append(value[end:match.start()], 'normal')
            token = match.group(0)
            if token[0] in _SUB:
                append(token.translate(_FROM_SUB), 'sub')
            elif token[0] in _SUP:
                append(token.translate(_FROM_SUP), 'sup')
            else:
                append(token[0], 'bar')
            end = match.end()
        append(value[end:], 'normal')

    value = str(text or '')
    end = 0
    for match in _EXPLICIT_SCRIPT.finditer(value):
        plain(value[end:match.start()])
        append(match[2].translate(_FROM_SUB).translate(_FROM_SUP),
               'sub' if match[1] == '_' else 'sup')
        end = match.end()
    plain(value[end:])
    return runs


def scientific_mathtext(text, *, bold=False, recognize_formulas=True):
    """Use upright Mathtext for scientific words, keeping prose in text mode."""
    # Preserve explicitly supplied Mathtext in custom Matplotlib titles.
    parts = re.split(r'((?<!\\)\$[^$\n]+(?<!\\)\$)', str(text or ''))
    if len(parts) > 1:
        return ''.join(part if index % 2 else scientific_mathtext(
            part, bold=bold, recognize_formulas=recognize_formulas)
                       for index, part in enumerate(parts))
    escapes = {'\\': r'\backslash ', '{': r'\{', '}': r'\}',
               '_': r'\_', '$': r'\$', '%': r'\%', ' ': r'\ ',
               '^': r'\text{^}'}

    def word(match):
        value = match.group(0)
        runs = scientific_runs(value, recognize_formulas=recognize_formulas)
        if not any(style != 'normal' for _, style in runs):
            return value.replace('$', r'\$')
        rendered = []
        for value, style in runs:
            part = ''.join(escapes.get(char, char) for char in value)
            if style == 'sub':
                rendered.append('_{' + part + '}')
            elif style == 'sup':
                rendered.append('^{' + part + '}')
            elif style == 'bar':
                rendered.append(r'\overline{' + part + '}')
            else:
                rendered.append(part)
        font = 'mathbf' if bold else 'mathrm'
        return '$\\' + font + '{' + ''.join(rendered) + '}$'

    return _LABEL_WORD.sub(word, str(text or ''))


class ScientificDraw:
    """Pillow drawing with matching measurement and baseline-aligned scripts.

    Non-text drawing and ordinary labels retain Pillow's native behavior.
    Styled runs preserve kerning within each run instead of drawing letters
    individually. All text coordinates use the same measured font metrics.
    """
    def __init__(self, image):
        self.draw = ImageDraw.Draw(image)

    def __getattr__(self, name):
        return getattr(self.draw, name)

    def _layout(self, text, font, anchor=None):
        runs = scientific_runs(text)
        if not any(style != 'normal' for _, style in runs):
            return None
        small = font.font_variant(size=max(7, round(font.size * 0.70)))
        ascent, descent = font.getmetrics()
        cursor = 0.0
        script_start = None
        placed, boxes = [], []
        for value, style in runs:
            if style in ('sub', 'sup'):
                if script_start is None:
                    script_start = cursor
                start = script_start
            else:
                script_start = None
                start = cursor
            face = small if style in ('sub', 'sup') else font
            shift = (font.size * 0.20 if style == 'sub' else
                     -font.size * 0.40 if style == 'sup' else 0)
            baseline = ascent + shift
            left, top, right, bottom = face.getbbox(value, anchor='ls')
            box = (start + left, baseline + top,
                   start + right, baseline + bottom)
            if style == 'bar':
                box = (box[0], box[1] - max(2, font.size * 0.10), box[2], box[3])
            boxes.append(box)
            placed.append((start, baseline, value, face, style, box))
            cursor = max(cursor, start + face.getlength(value))
        bounds = (min(b[0] for b in boxes), min(b[1] for b in boxes),
                  max(b[2] for b in boxes), max(b[3] for b in boxes))
        anchor = anchor or 'la'
        dx = {'l': 0, 'm': -cursor / 2, 'r': -cursor}[anchor[0]]
        dy = {'a': 0, 't': -bounds[1], 'm': -(ascent + descent) / 2,
              's': -ascent, 'b': -bounds[3], 'd': -(ascent + descent)}[anchor[1]]
        return placed, (bounds[0] + dx, bounds[1] + dy,
                        bounds[2] + dx, bounds[3] + dy), dx, dy, cursor

    def textbbox(self, xy, text, font=None, anchor=None, **kwargs):
        layout = self._layout(text, font, anchor)
        if layout is None:
            return self.draw.textbbox(xy, str(text), font=font, anchor=anchor, **kwargs)
        _, box, _, _, _ = layout
        return (math.floor(xy[0] + box[0]), math.floor(xy[1] + box[1]),
                math.ceil(xy[0] + box[2]), math.ceil(xy[1] + box[3]))

    def textlength(self, text, font=None, **kwargs):
        layout = self._layout(text, font)
        return (layout[4] if layout else self.draw.textlength(str(text), font=font, **kwargs))

    def text(self, xy, text, fill=None, font=None, anchor=None, **kwargs):
        layout = self._layout(text, font, anchor)
        if layout is None:
            return self.draw.text(xy, str(text), fill=fill, font=font, anchor=anchor, **kwargs)
        placed, _, dx, dy, _ = layout
        for x, baseline, value, face, style, box in placed:
            self.draw.text((xy[0] + dx + x, xy[1] + dy + baseline),
                           value, fill=fill, font=face, anchor='ls', **kwargs)
            if style == 'bar':
                self.draw.line((xy[0] + dx + box[0], xy[1] + dy + box[1],
                                xy[0] + dx + box[2], xy[1] + dy + box[1]),
                               fill=fill, width=max(1, round(font.size / 18)))


def scientific_text_image(text, font, fill=(0, 0, 0, 255), padding=6):
    """Tightly measure a label before rotating it, including descenders/scripts."""
    probe = ScientificDraw(Image.new('RGBA', (1, 1)))
    left, top, right, bottom = probe.textbbox((0, 0), text, font=font)
    image = Image.new('RGBA', (max(1, right - left + 2 * padding),
                               max(1, bottom - top + 2 * padding)))
    ScientificDraw(image).text((padding - left, padding - top), text, font=font, fill=fill)
    return image

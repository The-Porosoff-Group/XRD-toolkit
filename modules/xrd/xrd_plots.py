"""Publication-quality XRD refinement and candidate-preview figures."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter, MaxNLocator
import numpy as np

from modules.plot_style import arial_plot, scientific_mathtext
from .presentation import (
    clean_descriptive_text,
    enrich_phase_results,
    format_chemical_formula,
)


# Existing UI palette, retained for dark figures and candidate previews.
PHASE_COLORS = [
    '#f78166',  # coral red
    '#56d364',  # green
    '#e3b341',  # amber
    '#bc8cff',  # purple
    '#79c0ff',  # sky blue
    '#ffa657',  # orange
    '#ff7eb6',  # pink
    '#7ee787',  # mint
]

# Colorblind-accessible palette with sufficient contrast on a white page.
PUBLICATION_PHASE_COLORS = [
    '#0072b2',  # blue
    '#d55e00',  # vermilion
    '#009e73',  # bluish green
    '#cc79a7',  # reddish purple
    '#e69f00',  # orange
    '#56b4e9',  # sky blue
    '#7b61a8',  # violet
    '#8c6d31',  # brown
]

PLOT_THEMES = {
    'light': {
        'figure': '#ffffff', 'surface': '#ffffff', 'grid': '#d9dee7',
        'text': '#111827', 'muted': '#596273', 'observed': '#111827',
        'calculated': '#c1121f', 'residual': '#1d4ed8',
        'zero': '#4b5563', 'stats_face': '#f7f8fa',
        'phase_colors': PUBLICATION_PHASE_COLORS,
    },
    'dark': {
        'figure': '#0d1117', 'surface': '#161b22', 'grid': '#39414d',
        'text': '#f0f3f6', 'muted': '#a7b0bd', 'observed': '#e6edf3',
        'calculated': '#ff7b72', 'residual': '#79c0ff',
        'zero': '#8b949e', 'stats_face': '#1c2128',
        'phase_colors': PHASE_COLORS,
    },
}


def _phase_axis_label(label):
    """Return a single-line label for the phase-identification axis."""
    return str(label or '').strip()


def _mathtext_scientific_label(label, *, bold=False):
    """Typeset the explicit subscripts/bars from XRD presentation metadata."""
    return scientific_mathtext(label, bold=bold, recognize_formulas=False)


def _inclusive_two_theta_ticks(lower, upper):
    """Return whole-degree major ticks enclosing both data limits."""
    lower = float(lower)
    upper = float(upper)
    if upper <= lower:
        return np.asarray([np.floor(lower)])

    # Measured scans commonly begin or end a few floating-point increments
    # away from the requested range (for example, 10.0000000004 degrees).
    # Enclose the data in whole-degree display limits so those acquisition
    # artifacts never become long fractional labels on the exported plot.
    display_lower = np.floor(lower)
    display_upper = np.ceil(upper)
    if display_upper <= display_lower:
        display_upper = display_lower + 1

    locator = MaxNLocator(
        nbins=7, steps=[1, 2, 2.5, 5, 10], integer=True)
    epsilon = (display_upper - display_lower) * 1e-7
    interior = locator.tick_values(display_lower, display_upper)
    interior = interior[(interior > display_lower + epsilon)
                        & (interior < display_upper - epsilon)]
    interior = np.rint(interior[
        np.isclose(interior, np.rint(interior), atol=1e-9)])
    return np.unique(np.concatenate(
        ([display_lower], interior, [display_upper])))


@arial_plot
def make_xrd_plot(result, metadata, output_path, theme=None):
    """Render a fitted XRD pattern at its intended publication dimensions.

    Parameters
    ----------
    result : dict
        Completed refinement result.  Its numerical arrays are read-only here.
    metadata : dict
        Sample, wavelength, method, and optional ``plot_theme`` information.
        Set ``show_figure_title`` to False to omit the figure heading.
    output_path : str
        Destination PNG path.
    theme : {"light", "dark"}, optional
        Explicit figure theme. Light is the publication-export default.
    """
    theme = str(theme or metadata.get('plot_theme') or 'light').lower()
    if theme not in PLOT_THEMES:
        theme = 'light'
    palette = PLOT_THEMES[theme]

    # This enriches labels and tick metadata only.  Fit arrays and statistics
    # are not changed.
    enrich_phase_results(result)

    tt = np.asarray(result['tt'], dtype=float)
    y_obs = np.asarray(result['y_obs'], dtype=float)
    y_calc = np.asarray(result['y_calc'], dtype=float)
    y_bg = np.asarray(result['y_background'], dtype=float)
    resid = np.asarray(result['residuals'], dtype=float)
    phases = result.get('phase_results', []) or []
    phase_patterns = result.get('phase_patterns', []) or []
    n_phases = len(phases)
    stats = result['statistics']

    # Author at the final 6.5-inch manuscript width on a compact landscape
    # canvas so downstream software never has to stretch the export.
    extra_rows = max(n_phases - 2, 0)
    figure_height = min(4.95, 4.45 + 0.16 * extra_rows)
    # Each reflection row includes a heading above the rotated Miller indices.
    # Reserve physical space so adding phases cannot squeeze their text together.
    show_figure_title = metadata.get('show_figure_title', True)
    top_margin = 0.925 if show_figure_title else 0.975
    usable_fraction = (top_margin - 0.045) * 4 / (4 + 3 * 0.035)
    reflection_inches = 0.56 * max(n_phases, 1)
    figure_height = max(figure_height, (reflection_inches + 2.25) / usable_fraction)
    usable_inches = figure_height * usable_fraction
    tick_height = max(1.55, 6.55 * reflection_inches / (usable_inches - reflection_inches))

    fig = plt.figure(figsize=(6.5, figure_height),
                     facecolor=palette['figure'])
    grid = gridspec.GridSpec(
        4, 1, figure=fig, hspace=0.035,
        height_ratios=[4.8, 0.70, tick_height, 1.05])
    ax_main = fig.add_subplot(grid[0])
    ax_angle = fig.add_subplot(grid[1])
    ax_ticks = fig.add_subplot(grid[2], sharex=ax_main)
    ax_res = fig.add_subplot(grid[3], sharex=ax_main)

    surface = palette['surface']
    grid_color = palette['grid']
    text_color = palette['text']
    muted_color = palette['muted']
    phase_colors = palette['phase_colors']
    if theme == 'dark':
        observed_color = '#58a6ff'
        calculated_color = '#f0f6fc'
    else:
        observed_color = '#1d4ed8'
        calculated_color = '#111827'

    def style_axis(axis, show_xlabel=False, show_grid=True):
        axis.set_facecolor(surface)
        axis.tick_params(
            colors=text_color, labelsize=9.0, width=0.9, length=3.6,
            labelbottom=show_xlabel)
        axis.xaxis.label.set_color(text_color)
        axis.yaxis.label.set_color(text_color)
        for spine in axis.spines.values():
            spine.set_edgecolor(grid_color)
            spine.set_linewidth(0.9)
        if show_grid:
            axis.grid(True, color=grid_color, alpha=0.48, linewidth=0.65)
        else:
            axis.grid(False)
        axis.set_axisbelow(True)

    style_axis(ax_main, show_xlabel=True)
    style_axis(ax_ticks, show_grid=False)
    style_axis(ax_res)
    ax_main.tick_params(axis='x', pad=2)
    ax_angle.set_axis_off()
    ax_angle.text(
        0.5, 0.25, 'Diffraction angle, 2θ (degrees)',
        transform=ax_angle.transAxes, ha='center', va='center',
        fontsize=9.5, color=text_color)

    # Draw every phase against the fitted background.  These are true filled
    # component areas, not faint cumulative envelopes.
    for index, pattern in enumerate(phase_patterns[:n_phases]):
        color = phase_colors[index % len(phase_colors)]
        component = np.asarray(pattern, dtype=float)
        if component.size != y_bg.size:
            fitted_component = np.zeros_like(y_bg)
            copy_count = min(component.size, y_bg.size)
            fitted_component[:copy_count] = component[:copy_count]
            component = fitted_component
        phase_top = y_bg + np.maximum(component, 0)
        ax_main.fill_between(
            tt, y_bg, phase_top, color=color,
            alpha=0.42 if theme == 'dark' else 0.30,
            linewidth=0, zorder=1)
        ax_main.plot(tt, phase_top, color=color, linewidth=1.05,
                     alpha=0.98, zorder=2)

    ax_main.plot(tt, y_bg, color=muted_color, linewidth=0.95,
                 linestyle=(0, (4, 2)), alpha=0.95, zorder=2)
    ax_main.plot(tt, y_obs, color=observed_color, linewidth=1.05,
                 alpha=0.92, zorder=4)
    ax_main.plot(tt, y_calc, color=calculated_color, linewidth=1.45,
                 alpha=1.0, zorder=5)

    stats_text = (
        f"$R_{{\\mathrm{{wp}}}}$ {stats['Rwp']} %   "
        f"$R_{{\\mathrm{{p}}}}$ {stats['Rp']} %   "
        f"$\\chi^2$ {stats['chi2']}   GoF {stats['GoF']}"
    )
    custom_title = str(metadata.get('figure_title', '')).strip()
    sample_label = custom_title or (
        str(metadata.get('sample_id', 'Sample')).strip().replace('_', ' ')
        or 'Sample')
    title_fontsize = max(9.2, 12.4 - max(len(sample_label) - 42, 0) * 0.08)
    if show_figure_title:
        ax_main.set_title(
            scientific_mathtext(sample_label, bold=True),
            loc='center', pad=7, fontsize=title_fontsize,
            color=text_color, fontweight='bold')
    ax_main.text(
        0.995, 0.985, stats_text, transform=ax_main.transAxes,
        ha='right', va='top', fontsize=8.3, color=text_color,
        bbox=dict(boxstyle='round,pad=0.30', fc=palette['stats_face'],
                  ec=grid_color, alpha=0.90, linewidth=0.7), zorder=8)

    ax_main.set_ylabel('Intensity (a.u.)', fontsize=10.0,
                       color=text_color)
    ax_main.set_ylim(bottom=0)

    legend_handles = [
        Line2D([0], [0], color=observed_color, lw=1.45,
               label='Observed intensity'),
        Line2D([0], [0], color=calculated_color, lw=1.65,
               label='Calculated pattern'),
        Line2D([0], [0], color=muted_color, lw=1.2,
               ls=(0, (4, 2)), label='Fitted background'),
    ]
    for index, phase in enumerate(phases):
        color = phase_colors[index % len(phase_colors)]
        legend_handles.append(Patch(
            facecolor=color, edgecolor=color, alpha=0.70,
            label=_mathtext_scientific_label(phase['legend_label'])))
    figure_legend = ax_main.legend(
        handles=legend_handles, fontsize=8.2,
        ncol=1,
        facecolor=palette['stats_face'], edgecolor=grid_color,
        labelcolor=text_color, loc='upper right',
        bbox_to_anchor=(0.995, 0.875), frameon=True, fancybox=True,
        framealpha=0.88, borderpad=0.48, columnspacing=0.9,
        handlelength=1.65, handletextpad=0.48,
    )
    figure_legend.get_frame().set_linewidth(0.7)

    # Phase headings sit inside their rows, above the reflection labels.
    # No phase text extends into the left gutter.
    ax_ticks.set_ylim(0, max(n_phases, 1))
    ax_ticks.set_yticks([])
    for index, phase in enumerate(phases):
        color = phase_colors[index % len(phase_colors)]
        row_bottom = max(n_phases, 1) - index - 1
        ax_ticks.hlines(row_bottom + 0.10, tt.min(), tt.max(),
                        color=grid_color, linewidth=0.55, alpha=0.65)
        labeled_positions = set()
        for reflection in phase.get('tick_reflections', []) or []:
            position = float(reflection['two_theta'])
            label = reflection.get('label') or ''
            ax_ticks.vlines(position, row_bottom + 0.10, row_bottom + 0.32,
                            color=color, linewidth=1.45, alpha=1.0)
            if label:
                ax_ticks.text(
                    position, row_bottom + 0.36, _mathtext_scientific_label(label),
                    ha='center', va='bottom', rotation=60,
                    rotation_mode='anchor',
                    fontsize=7.0,
                    color=text_color, clip_on=True)
            labeled_positions.add(round(position, 3))
        for position in phase.get('tick_positions', []) or []:
            if round(float(position), 3) not in labeled_positions:
                ax_ticks.vlines(float(position), row_bottom + 0.10,
                                row_bottom + 0.32, color=color,
                                linewidth=1.45, alpha=1.0)
        phase_label = phase.get('tick_label') or clean_descriptive_text(
            phase.get('name', ''), fallback=f"Phase {index + 1}")
        ax_ticks.text(
            0.5, row_bottom + 0.96, _mathtext_scientific_label(
                _phase_axis_label(phase_label), bold=True),
            transform=ax_ticks.get_yaxis_transform(), ha='center', va='top',
            fontsize=7.5, fontweight='bold', linespacing=0.9,
            color=color, clip_on=True)

    ax_res.plot(tt, resid, color=palette['residual'], linewidth=0.95,
                alpha=1.0)
    ax_res.axhline(0, color=palette['zero'], linewidth=0.9,
                   linestyle=(0, (4, 2)), alpha=0.9)
    ax_res.fill_between(tt, resid, 0, where=(resid > 0),
                        color=palette['residual'], alpha=0.20)
    ax_res.fill_between(tt, resid, 0, where=(resid < 0),
                        color=phase_colors[0], alpha=0.20)
    ax_res.set_ylabel('Difference', fontsize=9.0, color=text_color)
    ax_res.tick_params(axis='y', labelsize=8.0)
    for axis in (ax_main, ax_res):
        axis.yaxis.set_label_coords(-0.085, 0.5)
        axis.yaxis.label.set_verticalalignment('center')

    two_theta_ticks = _inclusive_two_theta_ticks(tt.min(), tt.max())
    ax_main.set_xlim(two_theta_ticks[0], two_theta_ticks[-1])
    ax_main.set_xticks(two_theta_ticks)
    ax_main.xaxis.set_major_formatter(
        FuncFormatter(lambda value, _position: f'{value:.0f}'))
    fig.subplots_adjust(left=0.11, right=0.985, bottom=0.045,
                        top=top_margin)

    fig.savefig(
        output_path, dpi=300, facecolor=palette['figure'], edgecolor='none',
        metadata={
            'Title': f"{sample_label} XRD refinement",
            'Description': (
                f"{theme} theme; phase contributions shown as filled areas"
            ),
            'Software': 'Catalysis Data Toolkit',
        },
    )
    plt.close(fig)
    return output_path


@arial_plot
def make_candidate_preview(tt, y_obs, candidates, wavelength, output_path):
    """Render the pre-refinement candidate stick-pattern preview."""
    fig, ax = plt.subplots(figsize=(12, 5), facecolor='#0d1117')
    ax.set_facecolor('#161b22')
    ax.tick_params(colors='#e6edf3', labelsize=8)
    ax.set_xlabel('Diffraction angle, 2θ (degrees)', fontsize=9,
                  color='#e6edf3')
    ax.set_ylabel('Intensity (a.u.)', fontsize=9, color='#e6edf3')
    for spine in ax.spines.values():
        spine.set_edgecolor('#2d333b')
    ax.grid(True, color='#2d333b', alpha=0.4, linewidth=0.5)

    ax.plot(tt, y_obs, color='#58a6ff', linewidth=0.8, alpha=0.9,
            label='Observed intensity', zorder=3)

    for index, candidate in enumerate(candidates[:6]):
        color = PHASE_COLORS[index % len(PHASE_COLORS)]
        for stick in candidate.get('stick_pattern', []):
            ax.axvline(stick['two_theta'], color=color, linewidth=0.8,
                       alpha=0.5, linestyle='--', ymin=0, ymax=0.15)
        formula = _mathtext_scientific_label(
            format_chemical_formula(candidate.get('formula', '?')))
        ax.text(
            0.01 + index * 0.16, 0.97, f"●  {formula}",
            transform=ax.transAxes, ha='left', va='top', fontsize=7,
            color=color,
            bbox=dict(boxstyle='round,pad=0.2', fc='#1c2128',
                      ec=color, alpha=0.8))

    ax.set_title('Phase identification candidate overlay', color='#e6edf3',
                 fontsize=10, fontweight='bold')
    ax.set_xlim(np.min(tt), np.max(tt))
    plt.tight_layout()
    plt.savefig(output_path, dpi=130, bbox_inches='tight',
                facecolor='#0d1117')
    plt.close(fig)
    return output_path

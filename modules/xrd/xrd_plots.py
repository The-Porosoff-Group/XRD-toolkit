"""
modules/xrd/xrd_plots.py
Generates the XRD refinement output figure.

4-panel stacked layout:
  1. Raw data + total fit (main panel, tall)
  2. Per-phase tick marks (one row per phase)
  3. Per-phase calculated patterns (color-coded, stacked)
  4. Residuals (Yobs - Ycalc)
"""

import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

# Colour palette — distinct, accessible, dark-theme friendly
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


def make_xrd_plot(result, metadata, output_path):
    """
    result   : dict from lebail.run_lebail()
    metadata : dict with keys: sample_id, wavelength_label, etc.
    output_path : str — where to save the PNG
    """
    tt      = np.array(result['tt'])
    y_obs   = np.array(result['y_obs'])
    y_calc  = np.array(result['y_calc'])
    y_bg    = np.array(result['y_background'])
    resid   = np.array(result['residuals'])
    phases  = result['phase_results']
    n_ph    = len(phases)
    stats   = result['statistics']

    # ── Figure layout ────────────────────────────────────────────────────────
    # Rows: [main, ticks, residuals]  heights: [5, 1*n_phases, 1.5]
    n_rows    = 3
    h_ratios  = [5, 0.5 * max(n_ph, 1), 1.5]

    fig = plt.figure(figsize=(13, 4 + h_ratios[0] + h_ratios[1] + h_ratios[2]),
                     facecolor='#0d1117')
    gs  = gridspec.GridSpec(n_rows, 1, figure=fig,
                             hspace=0.0,
                             height_ratios=h_ratios)

    ax_main  = fig.add_subplot(gs[0])
    ax_ticks = fig.add_subplot(gs[1], sharex=ax_main)
    ax_res   = fig.add_subplot(gs[2], sharex=ax_main)

    BG   = '#0d1117'
    SURF = '#161b22'
    GRID = '#2d333b'
    TEXT = '#e6edf3'
    MUT  = '#7d8590'

    def style(ax, show_xlabel=False):
        ax.set_facecolor(SURF)
        ax.tick_params(colors=TEXT, labelsize=8)
        ax.xaxis.label.set_color(TEXT)
        ax.yaxis.label.set_color(TEXT)
        for sp in ax.spines.values():
            sp.set_edgecolor(GRID)
        ax.grid(True, color=GRID, alpha=0.4, linewidth=0.5)
        if not show_xlabel:
            plt.setp(ax.get_xticklabels(), visible=False)

    style(ax_main)
    style(ax_ticks)
    style(ax_res, show_xlabel=True)

    # ── Main panel: data + fit + background ──────────────────────────────────
    ax_main.plot(tt, y_obs, color='#58a6ff', linewidth=0.8,
                 alpha=0.85, label='$I_{obs}$', zorder=3)
    ax_main.plot(tt, y_calc, color='#ffffff', linewidth=1.2,
                 alpha=0.95, label='$I_{calc}$', zorder=4)
    ax_main.plot(tt, y_bg, color=MUT, linewidth=0.8,
                 linestyle='--', alpha=0.7, label='Background', zorder=2)

    # Per-phase patterns on main panel (stacked shaded fills)
    # Stack phases on top of each other so they don't overlap —
    # each phase's fill starts where the previous one ended.
    cumulative = np.array(y_bg, dtype=float)
    for i, (ph, pat) in enumerate(zip(phases, result['phase_patterns'])):
        color = PHASE_COLORS[i % len(PHASE_COLORS)]
        pat_arr = np.array(pat)
        new_top = cumulative + pat_arr
        ax_main.fill_between(tt, cumulative, new_top,
                              color=color, alpha=0.25, zorder=1)
        ax_main.plot(tt, new_top, color=color,
                     linewidth=0.7, alpha=0.6, zorder=2)
        cumulative = new_top

    # Stats annotation
    stats_str = (f"$R_{{wp}}$ = {stats['Rwp']}%   "
                 f"$R_p$ = {stats['Rp']}%   "
                 f"$\\chi^2$ = {stats['chi2']}   "
                 f"GoF = {stats['GoF']}")
    ax_main.text(0.99, 0.97, stats_str,
                 transform=ax_main.transAxes,
                 ha='right', va='top', fontsize=8,
                 color=TEXT, family='monospace',
                 bbox=dict(boxstyle='round,pad=0.4', fc='#1c2128',
                            ec=GRID, alpha=0.9))

    # Title
    lam_label = metadata.get('wavelength_label')
    if not lam_label:
        lam_label = f"λ={result.get('wavelength', 1.54056):.4f} Å"
    method_label = metadata.get('method', 'Le Bail')
    ax_main.set_title(
        f"{metadata.get('sample_id','Sample')}   ·   "
        f"{lam_label}   ·   {method_label} refinement",
        color=TEXT, fontsize=11, fontweight='bold', pad=10)

    ax_main.set_ylabel('Intensity (a.u.)', fontsize=9)
    ax_main.set_ylim(bottom=0)

    # Legend
    handles = [
        Line2D([0],[0], color='#58a6ff', lw=1.5, label='$I_{obs}$'),
        Line2D([0],[0], color='#ffffff', lw=1.5, label='$I_{calc}$'),
        Line2D([0],[0], color=MUT, lw=1, ls='--', label='Background'),
    ]
    for i, ph in enumerate(phases):
        c  = PHASE_COLORS[i % len(PHASE_COLORS)]
        wt = ph.get('weight_fraction_%', '')
        wt_err = ph.get('weight_fraction_err_%')
        if wt != '':
            wt_str = f"  {wt}"
            if wt_err is not None:
                wt_str += f" ± {wt_err}"
            wt_str += " wt%"
        else:
            wt_str = ''
        # ph['name'] already includes the space group for disambiguation
        label  = f"{ph['name']}{wt_str}"
        handles.append(Line2D([0],[0], color=c, lw=2, label=label))
    ax_main.legend(handles=handles, fontsize=7, ncol=min(len(handles), 4),
                   facecolor='#1c2128', edgecolor=GRID, labelcolor=TEXT,
                   loc='upper left')

    # ── Tick mark panel ──────────────────────────────────────────────────────
    ax_ticks.set_facecolor(SURF)
    ax_ticks.set_yticks([])
    ax_ticks.yaxis.set_visible(False)
    for sp in ax_ticks.spines.values():
        sp.set_edgecolor(GRID)

    y_positions = np.linspace(0.85, 0.15, max(n_ph, 1))
    for i, ph in enumerate(phases):
        color = PHASE_COLORS[i % len(PHASE_COLORS)]
        y_pos = y_positions[i]
        ticks = ph.get('tick_positions', [])
        for tt_tick in ticks:
            ax_ticks.axvline(tt_tick, ymin=y_pos-0.08, ymax=y_pos+0.08,
                              color=color, linewidth=1.0, alpha=0.8)
        # Phase label: name (already includes SG) + wt% ± err
        wt   = ph.get('weight_fraction_%', '')
        wt_err = ph.get('weight_fraction_err_%')
        if wt != '':
            wt_str = f"  {wt}"
            if wt_err is not None:
                wt_str += f" ± {wt_err}"
            wt_str += " wt%"
        else:
            wt_str = ''
        label  = f"{ph['name']}{wt_str}"
        ax_ticks.text(0.005, y_pos, label,
                      transform=ax_ticks.transAxes,
                      ha='left', va='center', fontsize=7,
                      color=color, fontweight='bold')

    # ── Residuals panel ──────────────────────────────────────────────────────
    ax_res.plot(tt, resid, color='#7d8590', linewidth=0.7, alpha=0.9)
    ax_res.axhline(0, color='#39d353', linewidth=0.8, linestyle='--', alpha=0.7)
    ax_res.fill_between(tt, resid, 0,
                         where=(resid > 0), color='#58a6ff', alpha=0.15)
    ax_res.fill_between(tt, resid, 0,
                         where=(resid < 0), color='#f78166', alpha=0.15)
    ax_res.set_ylabel('$I_{obs} - I_{calc}$', fontsize=8, color=TEXT)
    ax_res.set_xlabel('2θ (degrees)', fontsize=9, color=TEXT)

    # Align x axis
    ax_main.set_xlim(tt.min(), tt.max())

    plt.savefig(output_path, dpi=150, bbox_inches='tight',
                facecolor=BG)
    plt.close()
    return output_path


def make_candidate_preview(tt, y_obs, candidates, wavelength, output_path):
    """
    Quick preview plot: raw data with stick patterns of candidate phases overlaid.
    Used during phase identification before refinement.
    """
    fig, ax = plt.subplots(figsize=(12, 5), facecolor='#0d1117')
    ax.set_facecolor('#161b22')
    ax.tick_params(colors='#e6edf3', labelsize=8)
    ax.set_xlabel('2θ (degrees)', fontsize=9, color='#e6edf3')
    ax.set_ylabel('Intensity (a.u.)', fontsize=9, color='#e6edf3')
    for sp in ax.spines.values(): sp.set_edgecolor('#2d333b')
    ax.grid(True, color='#2d333b', alpha=0.4, linewidth=0.5)

    ax.plot(tt, y_obs, color='#58a6ff', linewidth=0.8, alpha=0.9,
            label='Data', zorder=3)

    ymax = y_obs.max()
    for i, cand in enumerate(candidates[:6]):
        color = PHASE_COLORS[i % len(PHASE_COLORS)]
        sticks = cand.get('stick_pattern', [])
        for s in sticks:
            ax.axvline(s['two_theta'], color=color,
                        linewidth=0.8, alpha=0.5, linestyle='--', ymin=0, ymax=0.15)
        # Label
        ax.text(0.01 + i*0.16, 0.97,
                f"●  {cand.get('formula','?')} [{cand.get('cod_id','')}]",
                transform=ax.transAxes, ha='left', va='top',
                fontsize=7, color=color,
                bbox=dict(boxstyle='round,pad=0.2', fc='#1c2128', ec=color, alpha=0.8))

    ax.set_title('Phase Identification Preview — Candidate Overlay',
                  color='#e6edf3', fontsize=10, fontweight='bold')
    ax.set_xlim(tt.min(), tt.max())
    plt.tight_layout()
    plt.savefig(output_path, dpi=130, bbox_inches='tight', facecolor='#0d1117')
    plt.close()
    return output_path

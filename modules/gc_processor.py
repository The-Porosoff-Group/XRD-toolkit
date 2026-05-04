"""
modules/gc_processor.py
Core GC data processing engine.
Loaded by app.py — do not run directly.
"""

import os, re, zipfile, xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import yaml


# ─────────────────────────────────────────────────────────────────────────────
# LOAD REACTION CONFIG
# ─────────────────────────────────────────────────────────────────────────────

def load_reaction_config(config_path):
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def list_reaction_configs(config_dir):
    configs = []
    for fname in sorted(os.listdir(config_dir)):
        if fname.endswith('.yaml') and fname != 'custom_template.yaml':
            path = os.path.join(config_dir, fname)
            with open(path, 'r') as f:
                cfg = yaml.safe_load(f)
            configs.append({
                'file':        fname,
                'name':        cfg.get('name', fname),
                'description': cfg.get('description', ''),
                'reactant':    cfg.get('reactant', ''),
                'inlet_species': cfg.get('inlet_species', []),
            })
    return configs


# ─────────────────────────────────────────────────────────────────────────────
# XLSX PARSER (direct XML — avoids openpyxl styling bug)
# ─────────────────────────────────────────────────────────────────────────────

def col_to_idx(letters):
    r = 0
    for ch in letters.upper():
        r = r * 26 + (ord(ch) - ord('A') + 1)
    return r - 1

def parse_xlsx(filepath):
    with zipfile.ZipFile(filepath) as z:
        with z.open('xl/sharedStrings.xml') as f:
            tree = ET.parse(f)
        ns = {'x': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
        strings = [
            si.find('.//x:t', ns).text if si.find('.//x:t', ns) is not None else ''
            for si in tree.findall('x:si', ns)
        ]
        with z.open('xl/worksheets/sheet.xml') as f:
            stree = ET.parse(f)

    rows = stree.findall('.//{http://schemas.openxmlformats.org/spreadsheetml/2006/main}row')

    def cell_val(c):
        t = c.get('t', '')
        v = c.find('{http://schemas.openxmlformats.org/spreadsheetml/2006/main}v')
        if v is None: return None
        return strings[int(v.text)] if t == 's' else v.text

    def row_dict(row):
        d = {}
        for c in row.findall('{http://schemas.openxmlformats.org/spreadsheetml/2006/main}c'):
            ref = c.get('r', '')
            letters = ''.join(ch for ch in ref if ch.isalpha())
            val = cell_val(c)
            if val is not None:
                d[col_to_idx(letters)] = val
        return d

    r0 = row_dict(rows[0])
    sequence_name = r0.get(1, 'Unknown')
    r2 = row_dict(rows[2])
    r3 = row_dict(rows[3])

    species_cols = {}
    sorted_species = sorted(r2.items())
    for col_idx, mtype in r3.items():
        if mtype == 'Amount':
            best, best_dist = None, 999
            for sc_idx, sc_name in sorted_species:
                d = abs(col_idx - sc_idx)
                if d <= 2 and d < best_dist:
                    best, best_dist = sc_name, d
            if best:
                species_cols[col_idx] = best

    injections = []
    for row in rows[4:]:
        d = row_dict(row)
        if not d: continue
        label = d.get(0, '')
        if not label: continue
        amounts = {}
        for cidx, sp in species_cols.items():
            val = d.get(cidx)
            if val is not None:
                try: amounts[sp] = float(val)
                except ValueError: pass
        m = re.search(r'(\d+)\s*$', label)
        injections.append({
            'label':     label,
            'inj_num':   int(m.group(1)) if m else None,
            'is_bypass': label.lower().startswith('bypass'),
            'amounts':   amounts,
        })

    return {'sequence_name': sequence_name, 'injections': injections}


# ─────────────────────────────────────────────────────────────────────────────
# MOLAR FLOW CALCULATIONS
# ─────────────────────────────────────────────────────────────────────────────

def find_ch4_tcd_key(species_config):
    for header, cfg in species_config.items():
        if cfg['label'] == 'CH4_TCD': return header
    return None

def find_ch4_fid_key(species_config):
    for header, cfg in species_config.items():
        if cfg['label'] == 'CH4': return header
    return None

def find_ar_key(species_config, is_label=None):
    for header, cfg in species_config.items():
        if cfg['label'] == 'Ar': return header
    return None

def compute_flows(amounts, F_Ar_sccm, species_config, use_ch4_bridge):
    ar_key = find_ar_key(species_config)
    C_Ar = amounts.get(ar_key) if ar_key else None
    if not C_Ar or C_Ar == 0:
        return {}

    ch4_ratio = None
    if use_ch4_bridge:
        tcd_key = find_ch4_tcd_key(species_config)
        fid_key = find_ch4_fid_key(species_config)
        c_tcd = amounts.get(tcd_key) if tcd_key else None
        c_fid = amounts.get(fid_key) if fid_key else None
        if c_tcd and c_fid and c_tcd > 0 and c_fid > 0:
            ch4_ratio = c_tcd / c_fid

    flows = {}
    for sp_header, cfg in species_config.items():
        C_A = amounts.get(sp_header)
        if C_A is None or C_A == 0: continue
        label = cfg['label']
        if cfg['det'] == 'TCD':
            flows[label] = F_Ar_sccm * (C_A / C_Ar)
        else:
            if use_ch4_bridge and ch4_ratio is not None:
                flows[label] = F_Ar_sccm * ch4_ratio * (C_A / C_Ar)
            elif not use_ch4_bridge:
                flows[label] = F_Ar_sccm * (C_A / C_Ar)
    return flows

def build_flow_table(data, F_Ar_sccm, species_config):
    ch4_tcd_key = find_ch4_tcd_key(species_config)
    ch4_fid_key = find_ch4_fid_key(species_config)
    has_bridge  = any(
        inj['amounts'].get(ch4_tcd_key) and inj['amounts'].get(ch4_fid_key)
        for inj in data['injections']
    )
    records = []
    for inj in data['injections']:
        flows = compute_flows(inj['amounts'], F_Ar_sccm, species_config, has_bridge)
        row = {'label': inj['label'], 'inj_num': inj['inj_num'],
               'is_bypass': inj['is_bypass']}
        row.update(flows)
        records.append(row)
    return pd.DataFrame(records), has_bridge

def get_cn(label, species_config):
    for cfg in species_config.values():
        if cfg['label'] == label: return cfg['cn']
    return 0

def calculate_results(df, reactant_label, F_reactant_inlet, species_config):
    df = df.copy()
    if reactant_label in df.columns:
        df['conversion'] = (F_reactant_inlet - df[reactant_label]) / F_reactant_inlet
    else:
        df['conversion'] = np.nan

    meta_cols = {'label', 'inj_num', 'is_bypass', 'conversion'}
    carbon_cols = [
        c for c in df.columns
        if c not in meta_cols
        and get_cn(c, species_config) > 0
        and c != reactant_label
    ]

    if carbon_cols:
        product_C  = sum(get_cn(c, species_config) * df[c].fillna(0) for c in carbon_cols)
        reactant_C = get_cn(reactant_label, species_config) * df[reactant_label].fillna(0) \
                     if reactant_label in df.columns else 0
        total_C_out = product_C + reactant_C
        with np.errstate(divide='ignore', invalid='ignore'):
            df_sel = pd.DataFrame({
                f'S_{c}': np.where(product_C > 0,
                    get_cn(c, species_config) * df[c].fillna(0) / product_C, np.nan)
                for c in carbon_cols
            }, index=df.index)
    else:
        total_C_out = pd.Series(np.nan, index=df.index)
        df_sel = pd.DataFrame(index=df.index)

    return df, df_sel, total_C_out, carbon_cols


# ─────────────────────────────────────────────────────────────────────────────
# PLOTTING
# ─────────────────────────────────────────────────────────────────────────────

def make_plots(df, df_sel, total_C_out, C_in_flow, reactant_label,
               ss_mask, metadata, carbon_cols, species_config, output_dir):

    fig = plt.figure(figsize=(15, 11))
    fig.patch.set_facecolor('#0f1117')
    fig.suptitle(
        f"{metadata['catalyst_id']}   ·   T={metadata.get('temperature','?')}   "
        f"P={metadata.get('pressure','?')}   GHSV={metadata.get('ghsv','?')}",
        fontsize=11, fontweight='bold', y=0.98, color='#e8e8e8')

    gs  = gridspec.GridSpec(2, 2, figure=fig, hspace=0.48, wspace=0.32)
    rxn = df[~df['is_bypass']].copy()
    inj = rxn['inj_num']

    ax_bg = '#161b22'
    ax_fg = '#e8e8e8'
    grid_c = '#2d333b'

    def style_ax(ax):
        ax.set_facecolor(ax_bg)
        ax.tick_params(colors=ax_fg, labelsize=8)
        ax.xaxis.label.set_color(ax_fg)
        ax.yaxis.label.set_color(ax_fg)
        ax.title.set_color(ax_fg)
        for spine in ax.spines.values():
            spine.set_edgecolor(grid_c)
        ax.grid(True, color=grid_c, alpha=0.6, linewidth=0.6)

    def shade(ax):
        if ss_mask[rxn.index].any():
            ss_inj = inj[ss_mask[rxn.index]]
            ax.axvspan(ss_inj.iloc[0] - 0.5, ss_inj.iloc[-1] + 0.5,
                       alpha=0.12, color='#39d353', zorder=0)

    colors = ['#58a6ff','#f78166','#56d364','#e3b341','#bc8cff',
              '#ff7b72','#79c0ff','#ffa657','#7ee787','#d2a8ff']

    # Panel 1: Molar flows
    ax1 = fig.add_subplot(gs[0, :])
    style_ax(ax1)
    plotted, ci = set(), 0
    priority = [reactant_label, 'CO2', 'CO', 'H2', 'CH4', 'C2H4', 'C2H6', 'C3H6', 'C3H8']
    for sp in priority + carbon_cols:
        if sp in rxn.columns and sp not in plotted and rxn[sp].notna().sum() > 0:
            ax1.plot(inj, rxn[sp], marker='o', markersize=3, linewidth=1.8,
                     label=sp, color=colors[ci % len(colors)])
            plotted.add(sp); ci += 1
    shade(ax1)
    ax1.set_xlabel('Injection number', fontsize=9)
    ax1.set_ylabel('Molar flow (sccm equiv.)', fontsize=9)
    ax1.set_title('Molar Flows vs. Injection Number')
    leg = ax1.legend(fontsize=7, ncol=6, loc='upper right',
                     facecolor='#1c2128', edgecolor=grid_c, labelcolor=ax_fg)

    # Panel 2: Conversion
    ax2 = fig.add_subplot(gs[1, 0])
    style_ax(ax2)
    if 'conversion' in rxn.columns and rxn['conversion'].notna().any():
        ax2.plot(inj, rxn['conversion'] * 100,
                 color='#58a6ff', marker='o', markersize=4, linewidth=2)
        shade(ax2)
        ss_rxn = ss_mask[rxn.index]
        if ss_rxn.any():
            ss_mean = rxn.loc[ss_rxn, 'conversion'].mean() * 100
            ax2.axhline(ss_mean, color='#39d353', linestyle='--', linewidth=1.5,
                        label=f'SS avg: {ss_mean:.2f}%')
            ax2.legend(fontsize=8, facecolor='#1c2128', edgecolor=grid_c, labelcolor=ax_fg)
    ax2.set_xlabel('Injection number', fontsize=9)
    ax2.set_ylabel('Conversion (%)', fontsize=9)
    ax2.set_title(f'{reactant_label} Conversion')

    # Panel 3: Selectivity bar
    ax3 = fig.add_subplot(gs[1, 1])
    style_ax(ax3)
    ss_rxn = ss_mask[rxn.index]
    if not df_sel.empty and ss_rxn.any():
        ss_sel = df_sel.loc[rxn.index][ss_rxn].mean() * 100
        ss_sel = ss_sel[ss_sel > 0.05].sort_values(ascending=False)
        if not ss_sel.empty:
            labels_bar = [s.replace('S_', '') for s in ss_sel.index]
            bar_colors = [colors[i % len(colors)] for i in range(len(ss_sel))]
            ax3.bar(range(len(ss_sel)), ss_sel.values, color=bar_colors, edgecolor='#0f1117', linewidth=0.5)
            ax3.set_xticks(range(len(ss_sel)))
            ax3.set_xticklabels(labels_bar, rotation=35, ha='right', fontsize=8)
            ax3.set_ylim(0, 115)
            if C_in_flow > 0:
                cb = (total_C_out[rxn.index][ss_rxn].mean() / C_in_flow) * 100
                cb_color = '#39d353' if 90 <= cb <= 110 else '#e3b341'
                ax3.text(0.98, 0.97, f'C balance: {cb:.1f}%',
                         transform=ax3.transAxes, ha='right', va='top', fontsize=9,
                         color=cb_color, fontweight='bold',
                         bbox=dict(boxstyle='round,pad=0.4', fc='#1c2128', ec=cb_color, alpha=0.9))
    ax3.set_ylabel('Selectivity (%)', fontsize=9)
    ax3.set_title('Carbon Selectivity (Steady State)')

    path = os.path.join(output_dir, 'gc_plots.png')
    plt.savefig(path, dpi=150, bbox_inches='tight', facecolor=fig.get_facecolor())
    plt.close()
    return path


# ─────────────────────────────────────────────────────────────────────────────
# SAVE OUTPUTS
# ─────────────────────────────────────────────────────────────────────────────

def save_outputs(df, df_sel, total_C_out, C_in_flow,
                 reactant_label, ss_mask, metadata, species_config, output_dir):
    rxn    = df[~df['is_bypass']]
    ss_rxn = ss_mask[rxn.index]
    row    = dict(metadata)

    row['n_bypass']   = int(df['is_bypass'].sum())
    row['n_reaction'] = int((~df['is_bypass']).sum())
    if ss_rxn.any():
        row['ss_inj_start'] = int(rxn.loc[ss_rxn, 'inj_num'].min())
        row['ss_inj_end']   = int(rxn.loc[ss_rxn, 'inj_num'].max())

    if 'conversion' in rxn.columns and ss_rxn.any():
        row['conversion_%']     = round(rxn.loc[ss_rxn, 'conversion'].mean() * 100, 2)
        row['conversion_std_%'] = round(rxn.loc[ss_rxn, 'conversion'].std()  * 100, 3)

    if not df_sel.empty and ss_rxn.any():
        for col, val in (df_sel.loc[rxn.index][ss_rxn].mean() * 100).items():
            sp = col.replace('S_', '')
            row[f'sel_{sp}_%'] = round(val, 2)

    if C_in_flow > 0 and ss_rxn.any():
        row['carbon_balance_%'] = round(
            (total_C_out[rxn.index][ss_rxn].mean() / C_in_flow) * 100, 2)

    summary = pd.DataFrame([row])
    summary_path = os.path.join(output_dir, 'gc_summary.csv')
    flows_path   = os.path.join(output_dir, 'gc_flows.csv')
    summary.to_csv(summary_path, index=False)
    df.to_csv(flows_path, index=False)
    return summary, summary_path, flows_path


# ─────────────────────────────────────────────────────────────────────────────
# MAIN ENTRY POINT (called by app.py)
# ─────────────────────────────────────────────────────────────────────────────

def run(filepath, output_dir, reaction_config, metadata, inlet_flows, ss_start, ss_end):
    """
    Main processing function called by the web app.
    Returns a dict with paths to output files and result summary.
    """
    os.makedirs(output_dir, exist_ok=True)

    species_config = reaction_config['species']
    reactant_label = reaction_config['reactant']

    # Find Ar inlet flow
    F_Ar = 0
    for sp, flow in inlet_flows.items():
        for header, cfg in species_config.items():
            if cfg['label'] == sp and cfg['label'] == 'Ar':
                F_Ar = flow
    if F_Ar == 0:
        F_Ar = inlet_flows.get('Ar', 15.0)

    # Carbon in
    C_in_flow = 0.0
    for sp, flow in inlet_flows.items():
        for cfg in species_config.values():
            if cfg['label'] == sp:
                C_in_flow += cfg['cn'] * flow
                break

    # Parse and compute
    data = parse_xlsx(filepath)
    df, has_bridge = build_flow_table(data, F_Ar, species_config)

    # Steady-state mask
    ss_mask = (~df['is_bypass']) & (df['inj_num'] >= ss_start) & (df['inj_num'] <= ss_end)

    df, df_sel, total_C_out, carbon_cols = calculate_results(
        df, reactant_label, inlet_flows.get(reactant_label, 0), species_config)

    plot_path = make_plots(
        df, df_sel, total_C_out, C_in_flow,
        reactant_label, ss_mask, metadata, carbon_cols, species_config, output_dir)

    summary, summary_path, flows_path = save_outputs(
        df, df_sel, total_C_out, C_in_flow,
        reactant_label, ss_mask, metadata, species_config, output_dir)

    # Build result dict for UI
    rxn    = df[~df['is_bypass']]
    ss_rxn = ss_mask[rxn.index]
    result = {
        'sequence_name':  data['sequence_name'],
        'n_bypass':       int(df['is_bypass'].sum()),
        'n_reaction':     int((~df['is_bypass']).sum()),
        'n_ss':           int(ss_rxn.sum()),
        'fid_bridge':     has_bridge,
        'plot_path':      plot_path,
        'summary_path':   summary_path,
        'flows_path':     flows_path,
        'output_dir':     output_dir,
    }
    if 'conversion' in rxn.columns and ss_rxn.any():
        result['conversion']     = round(rxn.loc[ss_rxn, 'conversion'].mean() * 100, 2)
        result['conversion_std'] = round(rxn.loc[ss_rxn, 'conversion'].std()  * 100, 3)
    if 'carbon_balance_%' in summary.columns:
        result['carbon_balance'] = round(summary['carbon_balance_%'].iloc[0], 2)

    sel_cols = [c for c in summary.columns if c.startswith('sel_') and summary[c].iloc[0] > 0.05]
    result['selectivities'] = {
        c.replace('sel_','').replace('_%',''): round(summary[c].iloc[0], 1)
        for c in sel_cols
    }

    return result

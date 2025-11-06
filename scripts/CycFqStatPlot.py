#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pypdf import PdfWriter, PdfReader
from matplotlib.ticker import MaxNLocator

# Global font settings
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['axes.labelweight'] = 'bold'

try:
    import plotly.graph_objects as go
    PLOTLY_AVAILABLE = True
except Exception:
    PLOTLY_AVAILABLE = False


def load_stats(csv_path: str) -> pd.DataFrame:
    # Support plain CSV and gzipped CSV (.csv.gz)
    df = pd.read_csv(csv_path, compression='infer')
    # Ensure expected columns exist
    expected = [
        'read_id', 'quality', 'length', 'has_adapter',
        'quality_pass', 'length_pass', 'final_output_count'
    ]
    missing = [c for c in expected if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in CSV: {missing}")
    # Normalize types
    df['has_adapter'] = df['has_adapter'].astype(str).str.lower().isin(['true', '1', 'yes'])
    df['quality_pass'] = df['quality_pass'].astype(str).str.lower().isin(['true', '1', 'yes'])
    df['length_pass'] = df['length_pass'].astype(str).str.lower().isin(['true', '1', 'yes'])
    df['final_output_count'] = pd.to_numeric(df['final_output_count'], errors='coerce').fillna(0).astype(int)
    df['quality'] = pd.to_numeric(df['quality'], errors='coerce')
    df['length'] = pd.to_numeric(df['length'], errors='coerce')
    return df


def plot_scatter(df: pd.DataFrame, min_quality: float, min_length: int, out_path: str):
    # Use original values; convert length to kb
    q = pd.to_numeric(df['quality'].replace([np.inf, -np.inf], np.nan), errors='coerce')
    l_kb = pd.to_numeric(df['length'].replace([np.inf, -np.inf], np.nan), errors='coerce') / 1000.0

    colors = np.where(df['has_adapter'], 'red', 'steelblue')

    # Layout with marginal density plots (top for X, right for Y)
    fig = plt.figure(figsize=(8, 6))
    gs = fig.add_gridspec(2, 2, width_ratios=[4, 1.2], height_ratios=[1.2, 4], wspace=0.05, hspace=0.05)
    ax = fig.add_subplot(gs[1, 0])
    ax_top = fig.add_subplot(gs[0, 0], sharex=ax)
    ax_right = fig.add_subplot(gs[1, 1], sharey=ax)

    # Main scatter
    ax.scatter(q, l_kb, c=colors, s=10, alpha=0.7, edgecolors='none')
    ax.set_xlabel('Quality (raw)', fontweight='bold')
    ax.set_ylabel('Length (kb)', fontweight='bold')

    # Fixed kb ticks (non-log)
    base_ticks = [1, 5, 10, 20, 25, 50, 100]
    try:
        y_max_data = np.nanmax(l_kb.values)
    except Exception:
        y_max_data = 0
    y_max_data = y_max_data if (np.isfinite(y_max_data) and y_max_data > 0) else 1
    ax.set_ylim(0, y_max_data)
    y_ticks = [t for t in base_ticks if t <= y_max_data]
    if y_ticks:
        ax.set_yticks(y_ticks)

    # Threshold lines
    ax.axvline(x=min_quality, color='orange', linestyle='--', linewidth=1.5, label=f'min_quality={min_quality}')
    ax.axhline(y=min_length/1000.0, color='green', linestyle='--', linewidth=1.5, label=f'min_length={min_length} bp')

    # Annotations for threshold lines
    try:
        x_max_data = np.nanmax(q.values)
    except Exception:
        x_max_data = min_quality if min_quality else 1.0
    y_max_data = ax.get_ylim()[1]
    min_length_kb = min_length / 1000.0

    ax.annotate(
        f"Quality ≥ {min_quality}",
        xy=(min_quality, y_max_data),
        xytext=(min_quality, y_max_data * 0.95),
        ha='center', va='top',
        color='orange',
        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.6, edgecolor='none'),
        arrowprops=dict(arrowstyle='-|>', color='orange', lw=1)
    )

    ax.annotate(
        f"Length ≥ {min_length_kb:g} kb",
        xy=(x_max_data, min_length_kb),
        xytext=(x_max_data * 0.98 if np.isfinite(x_max_data) else 0.98, min_length_kb),
        ha='right', va='center',
        color='green',
        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.6, edgecolor='none'),
        arrowprops=dict(arrowstyle='-|>', color='green', lw=1)
    )

    # Legend for adapter color at the top-left of the whole figure (place on ax_top)
    adapter_patch = plt.Line2D([0], [0], marker='o', color='w', label='Has adapter',
                               markerfacecolor='red', markersize=6)
    no_adapter_patch = plt.Line2D([0], [0], marker='o', color='w', label='No adapter',
                                  markerfacecolor='steelblue', markersize=6)
    ax_top.legend(handles=[adapter_patch, no_adapter_patch], loc='upper left', frameon=True)

    # --- Marginal density plots ---
    def kde_or_hist(xvals, grid_points=200):
        xvals = np.asarray(xvals)
        xvals = xvals[np.isfinite(xvals)]
        if xvals.size == 0:
            return np.array([0.0, 1.0]), np.array([0.0, 0.0])
        xmin, xmax = np.nanmin(xvals), np.nanmax(xvals)
        if not np.isfinite(xmin) or not np.isfinite(xmax) or xmax <= xmin:
            xmin, xmax = float(np.nanmin(xvals)), float(np.nanmax(xvals))
            if not np.isfinite(xmin) or not np.isfinite(xmax) or xmax <= xmin:
                xmax = xmin + 1.0
        grid = np.linspace(xmin, xmax, grid_points)
        n = xvals.size
        std = np.std(xvals, ddof=1) if n > 1 else 0.0
        if n > 10000 or std == 0.0:
            hist, edges = np.histogram(xvals, bins=min(50, max(10, int(np.sqrt(n)))), density=True)
            centers = 0.5 * (edges[:-1] + edges[1:])
            return centers, hist
        # Gaussian KDE with Silverman's rule
        bandwidth = 1.06 * std * (n ** (-1.0/5.0))
        diffs = (grid[:, None] - xvals[None, :]) / (bandwidth + 1e-12)
        dens = np.exp(-0.5 * diffs**2).sum(axis=1) / (n * (bandwidth + 1e-12) * np.sqrt(2*np.pi))
        return grid, dens

    # Top: density of quality (X)
    color_top = '#2c7fb8'   # blue
    gx, dx = kde_or_hist(q)
    ax_top.plot(gx, dx, color=color_top)
    ax_top.fill_between(gx, 0, dx, color=color_top, alpha=0.25)
    ax_top.set_yticks([])
    ax_top.set_ylabel('')
    ax_top.grid(False)
    for spine in ['right', 'top']:
        ax_top.spines[spine].set_visible(False)

    # Right: density of length (Y) plotted horizontally
    color_right = '#41ab5d'  # green
    gy, dy = kde_or_hist(l_kb)
    ax_right.plot(dy, gy, color=color_right)
    ax_right.fill_betweenx(gy, 0, dy, color=color_right, alpha=0.25)
    ax_right.set_xticks([])
    ax_right.set_xlabel('')
    ax_right.grid(False)
    for spine in ['top', 'right']:
        ax_right.spines[spine].set_visible(False)

    # Main axis cosmetics
    ax.grid(True, linestyle=':', alpha=0.4)
    fig.tight_layout()

    # Save PNG for convenience, and return figure for PDF assembly
    if out_path:
        fig.savefig(out_path, dpi=160)
    return fig


def plot_sankey(df: pd.DataFrame, out_path: str):
    if not PLOTLY_AVAILABLE:
        print('Plotly not available. Skipping Sankey plot.')
        return

    total = len(df)
    q_pass = int(df['quality_pass'].sum())
    q_fail = total - q_pass

    lp_mask = df['quality_pass']
    l_pass = int(df.loc[lp_mask, 'length_pass'].sum())
    l_fail = int(lp_mask.sum()) - l_pass

    # Length pass group
    len_pass_df = df.loc[lp_mask & df['length_pass']]
    a_found = int(len_pass_df['has_adapter'].sum())
    a_not = int(len_pass_df.shape[0]) - a_found

    # Adapter outcomes
    len_pass_adapter_df = len_pass_df.loc[len_pass_df['has_adapter']]
    op_af = int((len_pass_adapter_df['final_output_count'] > 0).sum())
    ff_af = int((len_pass_adapter_df['final_output_count'] == 0).sum())

    # No adapter => outputs produced (original read kept)
    op_anf = a_not  # all pass both, no adapter => kept

    # Layout: left Raw, middle Adapter check, right Clean/Filtered.
    # Aggregate quality/length failures into a single Raw→Filtered link.
    labels = [
        'Raw',
        'Adapter',
        'Clean', 'Filtered'
    ]
    idx = {label: i for i, label in enumerate(labels)}

    # Explicit positions
    x = [
        0.0,  # Raw
        0.5,  # Adapter (check stage)
        1.0,  # Clean
        1.0   # Filtered
    ]
    y = [
        0.5,  # Raw
        0.32, # Adapter (moved upward to avoid overlap)
        0.4,  # Clean
        0.6   # Filtered
    ]

    # values for flows
    agg_fail = q_fail + l_fail
    # show a thin line even when count is zero
    min_line_value = 1e-4
    agg_fail_display = agg_fail if agg_fail > 0 else min_line_value
    adapter_in = l_pass
    adapter_to_clean = op_af + op_anf
    adapter_to_filtered = ff_af

    sources = [
        idx['Raw'],            # raw -> filtered (quality+length fails aggregated)
        idx['Raw'],            # raw -> adapter
        idx['Adapter'],        # adapter -> clean
        idx['Adapter']         # adapter -> filtered
    ]
    targets = [
        idx['Filtered'],
        idx['Adapter'],
        idx['Clean'],
        idx['Filtered']
    ]
    # Display values (with thin line for zero) and raw values for percentages
    values_display = [
        agg_fail_display,
        adapter_in,
        adapter_to_clean,
        adapter_to_filtered
    ]
    values_raw = [
        agg_fail,
        adapter_in,
        adapter_to_clean,
        adapter_to_filtered
    ]
    values = values_display

    # Labels: show count + percentage of total
    def fmt_pct(v, tot):
        p = (v / tot) if tot > 0 else 0.0
        return f"{v} ({p*100:.1f}%)"
    link_labels = [fmt_pct(v, total) for v in values_raw]

    fig = go.Figure(data=[go.Sankey(
        arrangement="snap",
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=labels,
            x=x,
            y=y,
            color=[
                '#8da0cb', '#66c2a5', '#fc8d62',
                '#ffd92f', '#e78ac3',
                '#a6d854', '#b3b3b3',
                '#1b9e77', '#d95f02'
            ]
        ),
        link=dict(
            source=sources,
            target=targets,
            value=values,
            label=link_labels
        )
    )])
    # Add visible annotations near links (domain coordinates 0-1)
    ann_specs = [
        # (x, y, text); simplified to number + percentage to reduce clutter
        (0.50, 0.05, link_labels[0]),  # Raw -> Filtered (place at bottom)
        (0.25, 0.30, link_labels[1]),  # Raw -> Adapter (adjusted up)
        (0.78, 0.40, link_labels[2]),  # Adapter -> Clean (slightly up)
        (0.78, 0.68, link_labels[3])   # Adapter -> Filtered
    ]

    for x_pos, y_pos, text in ann_specs:
        fig.add_annotation(
            x=x_pos, y=y_pos,
            xref='paper', yref='paper',
            text=text,
            showarrow=False,
            font=dict(size=10, color='#222'),
            align='center',
            bgcolor='rgba(255,255,255,0.5)'
        )
    fig.update_layout(
        title_text="Filtering Flow: Raw → Adapter → Clean/Filtered",
        font=dict(family='Arial', size=12),
        width=800,
        height=600,
        margin=dict(l=10, r=10, t=40, b=10)
    )
    # Save HTML
    fig.write_html(out_path)
    return fig


def main():
    parser = argparse.ArgumentParser(description='Visualize read stats: scatter and Sankey.')
    parser.add_argument('--csv', required=True, help='Path to xxx.stat.csv')
    parser.add_argument('--min-quality', type=float, required=True, help='Minimum quality threshold used')
    parser.add_argument('--min-length', type=int, required=True, help='Minimum length threshold used')
    parser.add_argument('--out-dir', default=None, help='Output directory for plots')
    args = parser.parse_args()

    df = load_stats(args.csv)

    base = os.path.splitext(os.path.basename(args.csv))[0]
    out_dir = args.out_dir or os.path.dirname(os.path.abspath(args.csv))
    os.makedirs(out_dir, exist_ok=True)
    scatter_path = os.path.join(out_dir, f"{base}.scatter.png")
    sankey_path = os.path.join(out_dir, f"{base}.sankey.html")
    report_pdf_path = os.path.join(out_dir, f"{base}.report.pdf")

    scatter_fig = plot_scatter(df, args.min_quality, args.min_length, scatter_path)
    sankey_fig = plot_sankey(df, sankey_path) if PLOTLY_AVAILABLE else None

    # Export pages as separate PDFs (vector for Sankey) and then merge
    scatter_pdf_path = os.path.join(out_dir, f"{base}.scatter.page.pdf")
    scatter_fig.savefig(scatter_pdf_path)
    plt.close(scatter_fig)

    sankey_pdf_path = None
    if sankey_fig is not None:
        try:
            sankey_pdf_path = os.path.join(out_dir, f"{base}.sankey.page.pdf")
            sankey_fig.write_image(sankey_pdf_path, format='pdf')
        except Exception as e:
            print(f"Warning: Failed to export Sankey as vector PDF: {e}")

    # Merge into a single report
    writer = PdfWriter()
    # append scatter pages
    scatter_reader = PdfReader(scatter_pdf_path)
    for page in scatter_reader.pages:
        writer.add_page(page)
    # append sankey pages if available
    if sankey_pdf_path and os.path.exists(sankey_pdf_path):
        sankey_reader = PdfReader(sankey_pdf_path)
        for page in sankey_reader.pages:
            writer.add_page(page)
    with open(report_pdf_path, 'wb') as f_out:
        writer.write(f_out)

    print(f"Saved scatter: {scatter_path}")
    if PLOTLY_AVAILABLE:
        print(f"Saved sankey: {sankey_path}")
    else:
        print("Plotly not installed; sankey HTML not generated.")
    print(f"Saved report PDF: {report_pdf_path}")


if __name__ == '__main__':
    main()
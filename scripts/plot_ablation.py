#!/usr/bin/env python3
"""Plot ablation study results: precision and matching transcripts bar plots per parameter sweep."""

import argparse
import glob
import os
import re
import sys

import matplotlib.pyplot as plt
import numpy as np

# Default values matching main-ablation.nf
DEFAULTS = {
    'ir_part_v': 0.5,
    'ir_part_e': 0.5,
    'ir_full_v': 1.0,
    'ir_full_e': 0.5,
    'ir_full_i': 10.0,
    'umi_start_exon': 1,
    'umi_reads_bundle': 1,
    'umi_ratio_bundle': 0.0,
}

# Sorted by longest name first so greedy matching works
PARAM_NAMES = sorted(DEFAULTS.keys(), key=len, reverse=True)

# Most permissive (no-effect) values matching main-ablation.nf
# When a parameter is at its permissive value, it has no filtering effect.
PERMISSIVE = {
    'ir_part_v': 0,
    'ir_part_e': 0,
    'ir_full_v': 0,
    'ir_full_e': 0,
    'ir_full_i': 999,
    'umi_start_exon': 0,
    'umi_reads_bundle': 0,
    'umi_ratio_bundle': 0.0,
}


def parse_stats(filepath):
    """Parse a .compare.stats file.

    Returns (matching_count, precision) or None if parsing fails.
    """
    matching_count = 0
    precision = 0.0

    matching_key = 'Matching intron chains:'
    level_key = 'Intron chain level:'

    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith(matching_key):
                    matching_count = int(line.split(':')[1].strip())
                elif line.startswith(level_key):
                    parts = line.split('|')
                    if len(parts) >= 2:
                        precision = float(parts[1].strip())
    except (IOError, ValueError) as e:
        print(f"  Warning: failed to parse {filepath}: {e}", file=sys.stderr)
        return None

    return matching_count, precision


def split_label(dirname):
    """Split a directory label like 'ir_full_e_0.5' into (param_name, param_value).

    Uses known PARAM_NAMES to match the longest prefix.
    Returns (param_name, param_value) or None if no match.
    """
    for pname in PARAM_NAMES:
        prefix = pname + '_'
        if dirname.startswith(prefix):
            return pname, dirname[len(prefix):]
    return None


def discover_groups(results_dir, prefix='ablation'):
    """Glob ablation directories and group by param_name.

    Looks for directories matching {label}-{timestamp} under:
      {results_dir}/{prefix}-gffcpr/

    Returns dict: param_name -> list of (param_value, stats_dir)
    Also returns all_defaults dir and all_permissive dir separately.
    """
    gffcpr_base = os.path.join(results_dir, f'{prefix}-gffcpr')

    if not os.path.isdir(gffcpr_base):
        print(f"Error: gffcompare directory not found: {gffcpr_base}", file=sys.stderr)
        sys.exit(1)

    # Pattern: {label}-{YYYY-MM-DD_HH-MM-SS}
    timestamp_re = re.compile(r'^(.+)-(\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2})$')

    groups = {}
    all_defaults_dir = None
    all_permissive_dir = None
    for entry in sorted(os.listdir(gffcpr_base)):
        stats_dir = os.path.join(gffcpr_base, entry)
        if not os.path.isdir(stats_dir):
            continue

        m = timestamp_re.match(entry)
        if not m:
            print(f"  Warning: skipping unrecognized directory: {entry}", file=sys.stderr)
            continue

        label = m.group(1)

        # Skip baseline (old format)
        if label == 'baseline':
            continue

        if label == 'all_defaults':
            all_defaults_dir = stats_dir
            continue

        if label == 'all_permissive':
            all_permissive_dir = stats_dir
            continue

        parsed = split_label(label)
        if parsed is None:
            print(f"  Warning: could not match label to known param: {label}", file=sys.stderr)
            continue

        param_name, param_value = parsed

        if param_name not in groups:
            groups[param_name] = []
        groups[param_name].append((param_value, stats_dir))

    return groups, all_defaults_dir, all_permissive_dir


def _parse_stats_dir(stats_dir):
    """Parse all stats files in a directory. Returns list of (matching, precision)."""
    if stats_dir is None or not os.path.isdir(stats_dir):
        return []
    stats_files = glob.glob(os.path.join(stats_dir, '*.compare.stats'))
    results = []
    for sf in stats_files:
        result = parse_stats(sf)
        if result is not None:
            results.append(result)
    return results


def collect_data(groups, all_defaults_dir=None, all_permissive_dir=None):
    """For each (param_name, param_value), parse stats files and collect metrics.

    Returns:
        data: dict data[param_name][param_value] = [(matching, precision), ...]
        defaults_data: list of (matching, precision) or empty
        permissive_data: list of (matching, precision) or empty
    """
    data = {}
    defaults_data = []
    permissive_data = []

    for param_name, rows in groups.items():
        data[param_name] = {}
        for param_value, stats_dir in rows:
            cell_results = _parse_stats_dir(stats_dir)
            if not cell_results:
                print(f"  Warning: no stats for {param_name}={param_value}",
                      file=sys.stderr)
                continue
            data[param_name][param_value] = cell_results

    # Collect all_defaults data
    if all_defaults_dir is not None:
        defaults_data = _parse_stats_dir(all_defaults_dir)
        if defaults_data:
            print(f"  all_defaults: {len(defaults_data)} cells")

    # Collect all_permissive data
    if all_permissive_dir is not None:
        permissive_data = _parse_stats_dir(all_permissive_dir)
        if permissive_data:
            print(f"  all_permissive: {len(permissive_data)} cells")
    else:
        # Infer all_permissive from sweep data: any parameter at its permissive value
        for param_name in data:
            perm_cells = _find_permissive_cells(param_name, data[param_name])
            if perm_cells:
                permissive_data = perm_cells
                print(f"  all_permissive: inferred from {param_name} "
                      f"permissive value ({len(perm_cells)} cells)")
                break

    return data, defaults_data, permissive_data


def sort_values(param_name, values):
    """Sort parameter values numerically."""
    try:
        return sorted(values, key=lambda v: float(v))
    except ValueError:
        return sorted(values)


def is_default(param_name, param_value):
    """Check if param_value matches the default for param_name."""
    if param_name not in DEFAULTS:
        return False
    default = DEFAULTS[param_name]
    try:
        return float(param_value) == float(default)
    except ValueError:
        return str(param_value) == str(default)


def plot_boxplots(ax, values, all_cell_data, param_name, ylabel):
    """Draw boxplots with default value highlighted."""
    colors = []
    for v in values:
        if is_default(param_name, v):
            colors.append('#FF9500')
        else:
            colors.append('#468a71')

    bp = ax.boxplot(all_cell_data, positions=range(len(values)), patch_artist=True,
                    widths=0.6, zorder=2)
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    for median in bp['medians']:
        median.set_color('black')

    # Ensure y-axis largest tick exceeds max data value (including outliers)
    all_vals = [v for sublist in all_cell_data for v in sublist]
    if all_vals:
        data_max = max(all_vals)
        current_ylim = ax.get_ylim()
        if current_ylim[1] <= data_max:
            margin = (data_max - current_ylim[0]) * 0.05
            ax.set_ylim(current_ylim[0], data_max + margin)

    ax.set_xticks(range(len(values)))
    ax.set_xticklabels(values, fontsize=12, rotation=45, ha='right')
    ax.set_ylabel(ylabel, fontsize=14)
    ax.tick_params(axis='y', labelsize=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def _find_permissive_cells(param_name, param_data):
    """Find the cells for the most permissive (no-filtering) value of a parameter.

    Looks for the permissive value in the sweep data. If found, returns the cell
    results list; otherwise returns [].
    """
    if param_name not in PERMISSIVE:
        return []
    perm_val = PERMISSIVE[param_name]
    for v, cells in param_data.items():
        try:
            if float(v) == float(perm_val):
                return cells
        except ValueError:
            pass
    return []



def plot_defaults_vs_permissive(defaults_data, permissive_data, output_dir):
    """Generate boxplots comparing all_defaults vs all_permissive (no filtering)."""
    os.makedirs(output_dir, exist_ok=True)

    def_cells = defaults_data
    perm_cells = permissive_data

    if not def_cells and not perm_cells:
        return

    labels = []
    prec_data = []
    match_data = []
    colors = []

    if def_cells:
        labels.append('UMI-pruning')
        prec_data.append([c[1] for c in def_cells])
        match_data.append([c[0] for c in def_cells])
        colors.append('#FF9500')

    if perm_cells:
        labels.append('non-pruning')
        prec_data.append([c[1] for c in perm_cells])
        match_data.append([c[0] for c in perm_cells])
        colors.append('#468a71')

    if not labels:
        return

    # Precision plot
    fig, ax = plt.subplots(figsize=(4, 4.5))
    bp = ax.boxplot(prec_data, positions=range(len(labels)), patch_artist=True,
                    widths=0.5, zorder=2)
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    for median in bp['medians']:
        median.set_color('black')
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=12)
    ax.set_ylabel('precision (%)', fontsize=14)
    ax.set_ylim(0, 100)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'defaults_vs_permissive_precision.png'),
                dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'defaults_vs_permissive_precision.pdf'),
                dpi=300, bbox_inches='tight')
    plt.close(fig)

    # Matching count plot
    fig, ax = plt.subplots(figsize=(4, 4.5))
    bp = ax.boxplot(match_data, positions=range(len(labels)), patch_artist=True,
                    widths=0.5, zorder=2)
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    for median in bp['medians']:
        median.set_color('black')
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=12)
    ax.set_ylabel('# matching transcripts', fontsize=14)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'defaults_vs_permissive_matching.png'),
                dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'defaults_vs_permissive_matching.pdf'),
                dpi=300, bbox_inches='tight')
    plt.close(fig)

    # Console summary
    print(f"\ndefaults vs no filtering")
    print(f"  {'condition':<20} {'precision':>20} {'matching':>20} {'n_cells':>8}")
    for i, lbl in enumerate(labels):
        mp = np.mean(prec_data[i])
        sp = np.std(prec_data[i])
        mm = np.mean(match_data[i])
        sm = np.std(match_data[i])
        print(f"  {lbl:<20} {mp:>8.2f} +/- {sp:<8.2f}"
              f" {mm:>8.1f} +/- {sm:<8.1f} {len(prec_data[i]):>6}")

    print(f"  Saved to {output_dir}/defaults_vs_permissive_*.pdf")


def make_plots(data, output_dir, defaults_data=None, permissive_data=None):
    """Generate per-parameter box plots and summary figures."""
    os.makedirs(output_dir, exist_ok=True)
    if defaults_data is None:
        defaults_data = []
    if permissive_data is None:
        permissive_data = []

    def_cells = defaults_data

    for param_name in data:
        param_data = data[param_name]
        if not param_data:
            continue

        values = sort_values(param_name, list(param_data.keys()))
        if not values:
            continue

        # Collect per-value cell data
        all_precs = []
        all_matches = []
        n_cells_list = []

        for v in values:
            cells = param_data[v]
            precs = [c[1] for c in cells]
            matches = [c[0] for c in cells]
            all_precs.append(precs)
            all_matches.append(matches)
            n_cells_list.append(len(cells))

        # Plot A: Precision
        fig, ax = plt.subplots(figsize=(6, 4.5))
        plot_boxplots(ax, values, all_precs, param_name, 'precision (%)')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'{param_name}_precision.png'),
                    dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(output_dir, f'{param_name}_precision.pdf'),
                    dpi=300, bbox_inches='tight')
        plt.close(fig)

        # Plot B: Matching count
        fig, ax = plt.subplots(figsize=(6, 4.5))
        plot_boxplots(ax, values, all_matches, param_name, '# matching transcripts')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'{param_name}_matching.png'),
                    dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(output_dir, f'{param_name}_matching.pdf'),
                    dpi=300, bbox_inches='tight')
        plt.close(fig)

        # Console summary
        print(f"\n{param_name}")
        print(f"  {'value':<15} {'precision':>20} {'matching':>20} {'n_cells':>8}")
        for i, v in enumerate(values):
            default_marker = ' *' if is_default(param_name, v) else ''
            mp = np.mean(all_precs[i])
            sp = np.std(all_precs[i])
            mm = np.mean(all_matches[i])
            sm = np.std(all_matches[i])
            print(f"  {str(v):<15} {mp:>8.2f} +/- {sp:<8.2f}"
                  f" {mm:>8.1f} +/- {sm:<8.1f} {n_cells_list[i]:>6}{default_marker}")
        if def_cells:
            dp = np.mean([c[1] for c in def_cells])
            dm = np.mean([c[0] for c in def_cells])
            print(f"  {'all_defaults':<15} {dp:>8.2f}{'':>11}"
                  f" {dm:>8.1f}{'':>11} {len(def_cells):>6} (ref)")

    # Summary figure
    param_names = sorted(data.keys())
    if not param_names:
        return

    n_params = len(param_names)
    ncols = 2
    nrows = (n_params + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(14, 3.5 * nrows))
    axes = axes.flatten()

    for idx, param_name in enumerate(param_names):
        ax = axes[idx]
        param_data = data[param_name]
        if not param_data:
            ax.set_visible(False)
            continue

        values = sort_values(param_name, list(param_data.keys()))
        all_precs = [[c[1] for c in param_data[v]] for v in values]

        plot_boxplots(ax, values, all_precs, param_name, 'precision (%)')

    # Hide unused subplots
    for idx in range(n_params, len(axes)):
        axes[idx].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'ablation_summary.png'),
                dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'ablation_summary.pdf'),
                dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"\nSaved summary figure to {output_dir}/ablation_summary.pdf")


def main():
    parser = argparse.ArgumentParser(description='Plot ablation study results.')
    parser.add_argument('--results-dir', required=True,
                        help='Base results directory (e.g., results-hs-refseq-smart3)')
    parser.add_argument('--output-dir', default='./ablation_plots', help='Output directory for plots')
    parser.add_argument('--prefix', default='ablation',
                        help='Directory prefix (default: "ablation")')
    args = parser.parse_args()

    if not os.path.isdir(args.results_dir):
        print(f"Error: results directory not found: {args.results_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Discovering {args.prefix} directories in {args.results_dir} ...")
    groups, all_defaults_dir, all_permissive_dir = discover_groups(args.results_dir, prefix=args.prefix)
    print(f"Found {len(groups)} parameters: {', '.join(sorted(groups.keys()))}")
    for pname in sorted(groups.keys()):
        vals = [v for v, _ in groups[pname]]
        print(f"  {pname}: {vals}")
    if all_defaults_dir:
        print("  all_defaults: found")
    if all_permissive_dir:
        print("  all_permissive: found")

    print(f"\nCollecting stats ...")
    data, defaults_data, permissive_data = collect_data(
        groups, all_defaults_dir, all_permissive_dir)

    print(f"Generating plots in {args.output_dir} ...")
    make_plots(data, args.output_dir, defaults_data, permissive_data)

    # Defaults vs no-filtering comparison plot
    plot_defaults_vs_permissive(defaults_data, permissive_data, args.output_dir)

    print("\nDone.")


if __name__ == '__main__':
    main()

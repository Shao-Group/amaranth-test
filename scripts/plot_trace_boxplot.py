"""
Plot Nextflow execution trace: runtime and peak memory boxplots per assembler.

Usage:
  python plot_trace_boxplot.py --amaranth <trace.txt> --scallop2 <trace.txt> --stringtie <trace.txt> [--output-prefix trace_boxplot]
"""
import argparse
import matplotlib.pyplot as plt
import numpy as np
import re


def parse_time(s):
    """Parse Nextflow time string to seconds."""
    s = s.strip()
    match = re.match(r'([\d.]+)\s*(ms|s|m|h)', s)
    if not match:
        return None
    val, unit = float(match.group(1)), match.group(2)
    if unit == 'ms':
        return val / 1000
    if unit == 's':
        return val
    if unit == 'm':
        return val * 60
    if unit == 'h':
        return val * 3600
    return None


def parse_mem(s):
    """Parse Nextflow memory string to MB."""
    s = s.strip()
    match = re.match(r'([\d.]+)\s*(B|KB|MB|GB|TB)', s)
    if not match:
        return None
    val, unit = float(match.group(1)), match.group(2)
    scale = {'B': 1e-6, 'KB': 1e-3, 'MB': 1, 'GB': 1000, 'TB': 1e6}
    return val * scale[unit]


def extract_process_data(tracefile, process_prefix):
    """Extract realtime and peak_rss for rows matching process_prefix."""
    times, mems = [], []
    with open(tracefile) as f:
        header = f.readline().strip().split('\t')
        name_idx = header.index('name')
        rt_idx = header.index('realtime')
        rss_idx = header.index('peak_rss')
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) <= max(name_idx, rt_idx, rss_idx):
                continue
            if cols[name_idx].startswith(process_prefix):
                t = parse_time(cols[rt_idx])
                m = parse_mem(cols[rss_idx])
                if t is not None and m is not None:
                    times.append(t)
                    mems.append(m)
    return times, mems


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--amaranth',  required=True, help='Nextflow trace file for Amaranth run')
    parser.add_argument('--scallop2',  required=True, help='Nextflow trace file for Scallop2 run')
    parser.add_argument('--stringtie', required=True, help='Nextflow trace file for StringTie2 run')
    parser.add_argument('--output-prefix', default='trace_boxplot',
                        help='Output file prefix (default: trace_boxplot)')
    args = parser.parse_args()

    amr_time, amr_mem = extract_process_data(args.amaranth,  'RUNAMARANTH')
    sc2_time, sc2_mem = extract_process_data(args.scallop2,  'RUNSCALLOP2')
    stg_time, stg_mem = extract_process_data(args.stringtie, 'RUNSTRINGTIE')

    colors = ['#FF9500', '#468a71', '#a89c87']
    labels = ['Amaranth', 'Scallop2', 'StringTie2']

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7, 4))

    bp1 = ax1.boxplot([amr_time, sc2_time, stg_time], patch_artist=True, widths=0.6)
    for patch, color in zip(bp1['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    for element in ['whiskers', 'caps', 'medians']:
        for item in bp1[element]:
            item.set_color('black')
    ax1.set_xticklabels(labels, fontsize=11)
    ax1.set_ylabel('time (s)', fontsize=13)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.tick_params(axis='y', labelsize=11)

    bp2 = ax2.boxplot([amr_mem, sc2_mem, stg_mem], patch_artist=True, widths=0.6)
    for patch, color in zip(bp2['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    for element in ['whiskers', 'caps', 'medians']:
        for item in bp2[element]:
            item.set_color('black')
    ax2.set_xticklabels(labels, fontsize=11)
    ax2.set_ylabel('peak memory (MB)', fontsize=13)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.tick_params(axis='y', labelsize=11)

    plt.tight_layout()
    plt.savefig(f'{args.output_prefix}.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{args.output_prefix}.png', dpi=300, bbox_inches='tight')
    print(f'Saved {args.output_prefix}.pdf and {args.output_prefix}.png')

    for label, times, mems in zip(labels, [amr_time, sc2_time, stg_time], [amr_mem, sc2_mem, stg_mem]):
        print(f'\n{label}: n={len(times)}')
        print(f'  realtime: {np.median(times):.1f}s (median), {np.mean(times):.1f}s (mean)')
        print(f'  peak_rss: {np.median(mems):.1f} MB (median), {np.mean(mems):.1f} MB (mean)')

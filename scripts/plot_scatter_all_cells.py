"""
Scatter plot of per-cell assembly performance: # matching transcripts vs precision.

Reads TSV files produced by STATS2TSV (modules.nf) at intron-chain or transcript level,
and plots one point per cell for each assembler.

Usage:
  python plot_scatter_all_cells.py --recall <number.tsv> --precision <precision.tsv> [--amr-row <row_name>] [--output-prefix all.scatter]
"""
import argparse
import matplotlib.pyplot as plt
import numpy as np

COLORS = ['#FF9500', '#468a71', '#a89c87']  # Amaranth, Scallop2, StringTie2

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--recall',    required=True, help='TSV file with matching transcript counts (from STATS2TSV)')
    parser.add_argument('--precision', required=True, help='TSV file with precision values (from STATS2TSV)')
    parser.add_argument('--output-prefix', default='all.scatter',
                        help='Output file prefix (default: all.scatter)')
    args = parser.parse_args()

    amr = "amr"
    sc2 = "sc2"
    stg = "stg"

    amr_corr, amr_prec = [], []
    sc2_corr, sc2_prec = [], []
    stg_corr, stg_prec = [], []

    with open(args.recall, 'r') as file:
        for line in file:
            name = line.split("\t")[0]
            if name == amr:
                amr_corr = [int(x) for x in line.split("\t")[2:]]
            if name == sc2:
                sc2_corr = [int(x) for x in line.split("\t")[2:]]
            if name == stg:
                stg_corr = [int(x) for x in line.split("\t")[2:]]

    with open(args.precision, 'r') as file:
        for line in file:
            name = line.split("\t")[0]
            if name == amr:
                amr_prec = [float(x) for x in line.split("\t")[2:]]
            if name == sc2:
                sc2_prec = [float(x) for x in line.split("\t")[2:]]
            if name == stg:
                stg_prec = [float(x) for x in line.split("\t")[2:]]

    assert len(amr_corr) == len(sc2_corr) == len(stg_corr) == len(amr_prec) == len(sc2_prec) == len(stg_prec)

    plt.figure(figsize=(5, 5))
    correct_values_all = [amr_corr, sc2_corr, stg_corr]
    precision_values_all = [amr_prec, sc2_prec, stg_prec]
    methods = [amr, sc2, stg]
    method_names = {sc2: "Scallop2", stg: "StringTie2", amr: "Amaranth"}

    for idx, method in enumerate(methods):
        correct_values = correct_values_all[idx]
        precision_values = precision_values_all[idx]

        print(f"\n{'='*60}")
        print(f"Method: {method}")
        print(f"{'='*60}")
        print(f"Recall (# matching transcripts):")
        print(f"  Average: {np.mean(correct_values):.2f} ± {np.std(correct_values):.2f}")
        print(f"  Range: {min(correct_values)} - {max(correct_values)}")
        print(f"Precision (%):")
        print(f"  Average: {np.mean(precision_values):.2f} ± {np.std(precision_values):.2f}")
        print(f"  Range: {min(precision_values):.2f} - {max(precision_values):.2f}")

        methodname = method_names[method]
        plt.scatter(correct_values, precision_values, color=COLORS[idx], s=10, zorder=2, label=methodname, alpha=0.6)

    plt.xlabel('# matching transcripts', fontsize=14)
    plt.ylabel('precision (%)', fontsize=14)
    plt.ylim(0, 100)
    plt.xticks(fontsize=12, rotation=45, ha='right')
    plt.yticks(fontsize=12, rotation=90, va='center')
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()

    plt.savefig(f"{args.output_prefix}.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{args.output_prefix}.png", dpi=300, bbox_inches='tight')

    handles, labels = plt.gca().get_legend_handles_labels()
    fig_legend = plt.figure(figsize=(3, 0.5))
    fig_legend.legend(handles, labels, loc='center', ncol=3, fontsize=12, frameon=True, fancybox=True)
    fig_legend.savefig(f"{args.output_prefix}_legend.pdf", dpi=300, bbox_inches='tight')
    fig_legend.savefig(f"{args.output_prefix}_legend.png", dpi=300, bbox_inches='tight')
    plt.close(fig_legend)

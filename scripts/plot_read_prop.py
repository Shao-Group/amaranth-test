"""
Plot read proportion figures from RSeQC read_distribution and infer_experiment outputs.

Generates stacked bar plots and boxplots comparing internal reads vs UMI-tagged reads
by genomic feature (exonic, intronic, intergenic) and strand type (FR, RF, unknown).

Usage:
  python plot_read_prop.py --read-dist-dir <read_dist_dir> --strand-dir <strand_dir> [--output-dir <output_dir>]
"""
import argparse
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy


def parse_read_dist(filepath):
    result = {}
    with open(filepath) as f:
        lines = f.readlines()
    for line in lines:
        if line.strip().startswith("Total"):
            parts = line.rsplit(maxsplit=1)
            value = parts[-1]
            key = parts[0].strip()
            result[key] = int(value)
        elif line.strip().startswith("==============="):
            continue
        elif line.strip().startswith("Group"):
            continue
        else:
            parts = line.split()
            if len(parts) == 4:
                group = parts[0] + "-Tag_count"
                result[group] = int(parts[2])
    return result

def parse_libtype_file(filepath):
    result = {}
    with open(filepath) as f:
        lines = f.readlines()
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith("This is"):
            result["is_paired_end"] = "PairEnd" in line
        elif line.startswith("Fraction of reads failed to determine"):
            result["failed_fraction"] = float(line.split(": ")[1])
        elif line.startswith("Fraction of reads explained by \"1++,1--,2+-,2-+\""):
            result["fr_fraction"] = float(line.split(": ")[1])
        elif line.startswith("Fraction of reads explained by \"1+-,1-+,2++,2--\""):
            result["rf_fraction"] = float(line.split(": ")[1])
        elif line in ["FR", "RF", "Unkown"]:
            result["library_type"] = line
    return result


def flatten_total_counts(read_dist_all_files, use_tags=True):
    assert len(read_dist_all_files) % 2 == 0
    total_counts = [-1] * (len(read_dist_all_files) // 2)
    for i in range(len(read_dist_all_files)//2):
        id1, read_dist1 = read_dist_all_files[i * 2]
        id2, read_dist2 = read_dist_all_files[i * 2 + 1]
        id = id1.rsplit('.', 1)[0]
        assert id1 != id2 and id1 == id + ".internal" and id2 == id + ".umied"
        if use_tags:
            n1 = read_dist1["Total Tags"]
            n2 = read_dist2["Total Tags"]
        else:
            n1 = read_dist1["Total Reads"]
            n2 = read_dist2["Total Reads"]
        total_counts[i] = (n1, n2, n1 + n2)
    return total_counts

def plot_read_num_comp(total_counts, use_tags=False, output_dir='.'):
    total_counts = [[y / x[-1] for y in x] for x in total_counts]
    total_counts = sorted(total_counts, key=lambda x: x[0])
    comp1 = [x[0] for x in total_counts]
    comp2 = [x[1] for x in total_counts]

    fig = plt.figure(figsize=(5, 5))
    plt.tight_layout()
    plt.bar(range(len(comp1)), comp1, width=1.0, label='Internal reads', color='#FBC2B5', edgecolor='none', linewidth=0)
    plt.bar(range(len(comp1)), comp2, width=1.0, bottom=comp1, label='UMI reads', color='#95F9E3', edgecolor='none', linewidth=0)

    plt.margins(x=0)
    avg_comp1 = np.mean(comp1)
    plt.axhline(y=avg_comp1, color='black', linestyle='--', alpha=0.5)
    plt.text(len(comp1), avg_comp1, f'Mean: {int(avg_comp1*100)}%', verticalalignment='bottom', horizontalalignment='right', fontsize=12)

    plt.xlabel('Cells (n = {})'.format(len(total_counts)), fontsize=12)
    plt.xticks([], fontsize=12)
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1], ['0%', '20%', '40%', '60%', '80%', '100%'], fontsize=12, rotation=90, va='center')
    plt.ylim(0, 1)
    plt.tight_layout(pad=1.2)

    if use_tags:
        name_png = os.path.join(output_dir, "read_num_comp_classify_by_tags.png")
        name_pdf = os.path.join(output_dir, "read_num_comp_classify_by_tags.pdf")
        legend_name_png = os.path.join(output_dir, "read_num_comp_classify_by_tags_legend.png")
        legend_name_pdf = os.path.join(output_dir, "read_num_comp_classify_by_tags_legend.pdf")
    else:
        name_png = os.path.join(output_dir, "read_num_comp_classify_by_reads.png")
        name_pdf = os.path.join(output_dir, "read_num_comp_classify_by_reads.pdf")
        legend_name_png = os.path.join(output_dir, "read_num_comp_classify_by_reads_legend.png")
        legend_name_pdf = os.path.join(output_dir, "read_num_comp_classify_by_reads_legend.pdf")

    plt.savefig(name_png, dpi=600, bbox_inches='tight')
    plt.savefig(name_pdf, bbox_inches='tight')

    handles, labels = plt.gca().get_legend_handles_labels()
    handles = [handles[1], handles[0]]
    labels = [labels[1], labels[0]]
    fig_legend = plt.figure(figsize=(2, 0.5))
    fig_legend.legend(handles, labels, loc='center', ncol=2, fontsize=12, frameon=True, fancybox=True, handlelength=1.0)
    fig_legend.savefig(legend_name_png, dpi=600, bbox_inches='tight')
    fig_legend.savefig(legend_name_pdf, bbox_inches='tight')
    plt.close(fig_legend)
    plt.clf()

def plot_read_num_comp_error_bar(total_counts, output_dir='.'):
    total_counts = [[y / x[-1] for y in x] for x in total_counts]
    total_counts = sorted(total_counts, key=lambda x: x[0])
    comp1 = [x[0] for x in total_counts]
    comp2 = [x[1] for x in total_counts]

    avg_comp1 = sum(comp1) / len(comp1)
    avg_comp2 = sum(comp2) / len(comp2)
    assert avg_comp1 + avg_comp2 < 1.0 + 0.001 and avg_comp1 + avg_comp2 > 1.0 - 0.001

    fig = plt.figure(figsize=(1.3, 4))
    std_comp1 = np.std(comp1)

    plt.bar([0], [avg_comp1] * 1, label='Internal reads', color='#FBC2B5', zorder=1, edgecolor='none', linewidth=0)
    plt.bar([0], [avg_comp2] * 1, bottom=[avg_comp1] * 1, label='UMI reads', color='#95F9E3', zorder=1, edgecolor='none', linewidth=0)
    plt.errorbar([0], [avg_comp1] * 1, yerr=std_comp1, color='black', capsize=5, fmt='none', zorder=10)

    plt.xlabel('Cells', fontsize=12)
    plt.xticks([], fontsize=12)
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1], ['0%', '20%', '40%', '60%', '80%', '100%'], fontsize=12, rotation=90, va='center')
    plt.tight_layout(pad=1.2)

    plt.savefig(os.path.join(output_dir, "read_num_comp_error_bar.png"), dpi=600, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, "read_num_comp_error_bar.pdf"), bbox_inches='tight')

    handles, labels = plt.gca().get_legend_handles_labels()
    handles = [handles[1], handles[0]]
    labels = [labels[1], labels[0]]
    fig_legend = plt.figure(figsize=(2, 0.5))
    fig_legend.legend(handles, labels, loc='center', ncol=2, fontsize=12, frameon=True, fancybox=True, handlelength=1.0)
    fig_legend.savefig(os.path.join(output_dir, "read_num_comp_error_bar_legend.png"), dpi=600, bbox_inches='tight')
    fig_legend.savefig(os.path.join(output_dir, "read_num_comp_error_bar_legend.pdf"), bbox_inches='tight')
    plt.close(fig_legend)
    plt.clf()

def flatten_read_distributions(read_dist_all_files):
    assert len(read_dist_all_files) % 2 == 0
    total_counts = [-1] * (len(read_dist_all_files) // 2)
    for i in range(len(read_dist_all_files)//2):
        id1, read_dist1 = read_dist_all_files[i * 2]
        id2, read_dist2 = read_dist_all_files[i * 2 + 1]
        id = id1.rsplit('.', 1)[0]
        assert id1 != id2 and id1 == id + ".internal" and id2 == id + ".umied"
        exon1 = sum([v for k,v in read_dist1.items() if "Exons" in k])
        intron1 = sum([v for k,v in read_dist1.items() if "Introns" in k])
        intergenic1 = sum([v for k,v in read_dist1.items() if "TSS_up_10kb" in k or "TES_down_10kb" in k])
        unassigned1 = read_dist1["Total Tags"] - read_dist1["Total Assigned Tags"]
        exon2 = sum([v for k,v in read_dist2.items() if "Exons" in k])
        intron2 = sum([v for k,v in read_dist2.items() if "Introns" in k])
        intergenic2 = sum([v for k,v in read_dist2.items() if "TSS_up_10kb" in k or "TES_down_10kb" in k])
        unassigned2 = read_dist2["Total Tags"] - read_dist2["Total Assigned Tags"]
        assert exon1 + intron1 + intergenic1 == read_dist1["Total Assigned Tags"]
        assert exon2 + intron2 + intergenic2 == read_dist2["Total Assigned Tags"]
        all = exon1 + intron1 + intergenic1 + unassigned1 + exon2 + intron2 + intergenic2 + unassigned2
        total_counts[i] = (exon1, intron1, intergenic1, unassigned1, exon2, intron2, intergenic2, unassigned2, all)
    return total_counts

def flatten_strand_counts(strand_dist_all_files):
    assert len(strand_dist_all_files) % 2 == 0
    total_counts = [-1] * (len(strand_dist_all_files) // 2)
    for i in range(len(strand_dist_all_files)//2):
        id1, read_dist1 = strand_dist_all_files[i * 2]
        id2, read_dist2 = strand_dist_all_files[i * 2 + 1]
        id = id1.rsplit('.', 2)[0]
        assert id1 != id2 and id1 == id + ".internal.libtype" and id2 == id + ".umied.libtype"
        fr1 = read_dist1["fr_fraction"]
        rf1 = read_dist1["rf_fraction"]
        uk1 = read_dist1["failed_fraction"]
        fr2 = read_dist2["fr_fraction"]
        rf2 = read_dist2["rf_fraction"]
        uk2 = read_dist2["failed_fraction"]
        total_counts[i] = (fr1, rf1, uk1, fr2, rf2, uk2)
    return total_counts

def plot_read_dist_types(total_counts, output_dir='.'):
    total_counts = [[y / x[-1] for y in x] for x in total_counts]
    total_counts = sorted(total_counts, key=lambda x: sum(x[0:4]))
    fig = plt.figure(figsize=(5, 5))
    plt.tight_layout()
    colors = ['#FBC2B5', '#FFA8A9', '#F786AA', '#A14A76',
                '#95F9E3', '#69EBD0', '#49D49D', '#558564']
    labels = ['Internal Exonic', 'Internal Intronic', 'Internal Intergenic', 'Internal Unassigned',
              'UMI Exonic', 'UMI Intronic', 'UMI Intergenic', 'UMI Unassigned']

    bottom = np.zeros(len(total_counts))
    for i in range(8):
        values = [x[i] for x in total_counts]
        plt.bar(range(len(total_counts)), values, width=1.0, bottom=bottom, label=labels[i], color=colors[i], edgecolor='none', linewidth=0)
        bottom += values

    plt.margins(x=0)
    plt.xlabel('Cells (n = {})'.format(len(total_counts)), fontsize=12)
    plt.ylim(0, 1)
    plt.xticks([], fontsize=12)
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1], ['0%', '20%', '40%', '60%', '80%', '100%'], fontsize=12, rotation=90, va='center')
    plt.tight_layout(pad=1.2)

    plt.savefig(os.path.join(output_dir, "read_dist_types.png"), dpi=600, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, "read_dist_types.pdf"), bbox_inches='tight')

    handles, labels = plt.gca().get_legend_handles_labels()
    reordered_handles = []
    reordered_labels = []
    for i in range(4):
        reordered_handles.append(handles[i])
        reordered_handles.append(handles[i+4])
        reordered_labels.append(labels[i])
        reordered_labels.append(labels[i+4])
    fig_legend = plt.figure(figsize=(12, 1.2))
    fig_legend.legend(reordered_handles, reordered_labels, loc='center', ncol=4, fontsize=12, frameon=True, fancybox=True, handlelength=1.0)
    fig_legend.savefig(os.path.join(output_dir, "read_dist_types_legend.png"), dpi=600, bbox_inches='tight')
    fig_legend.savefig(os.path.join(output_dir, "read_dist_types_legend.pdf"), bbox_inches='tight')
    plt.close(fig_legend)
    plt.clf()

def get_strand_dist(original_total_counts, strand_counts):
    assert len(original_total_counts) == len(strand_counts)
    total_counts = []
    for total, strand in zip(original_total_counts, strand_counts):
        internal_fr, internal_rf, internal_unkown = total[0] * strand[0], total[0] * strand[1], total[0] * strand[2]
        umi_fr, umi_rf, umi_unkown = total[1] * strand[3], total[1] * strand[4], total[1] * strand[5]
        sum = internal_fr + internal_rf + internal_unkown + umi_fr + umi_rf + umi_unkown
        assert sum / total[-1] < 1.0 + 0.001 and sum / total[-1] > 1.0 - 0.001
        total_counts.append((internal_fr/sum, internal_rf/sum, internal_unkown/sum, umi_fr/sum, umi_rf/sum, umi_unkown/sum))
    return total_counts

def plot_strand_dist(total_counts, output_dir='.'):
    total_counts = sorted(total_counts, key=lambda x: sum(x[0:3]))
    fig = plt.figure(figsize=(5, 5))
    plt.tight_layout()
    colors = ['#FBC2B5', '#FFA8A9', '#A14A76',
                '#95F9E3', '#69EBD0', '#558564']
    labels = ['Internal FR', 'Internal RF', 'Internal Unknown',
              'UMI FR', 'UMI RF', 'UMI Unknown']

    bottom = np.zeros(len(total_counts))
    for i in range(6):
        values = [x[i] for x in total_counts]
        plt.bar(range(len(total_counts)), values, width=1.0, bottom=bottom, label=labels[i], color=colors[i], edgecolor='none', linewidth=0)
        bottom += values

    plt.margins(x=0)
    plt.xlabel('Cells (n = {})'.format(len(total_counts)), fontsize=12)
    plt.ylim(0, 1)
    plt.xticks([], fontsize=12)
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1], ['0%', '20%', '40%', '60%', '80%', '100%'], fontsize=12, rotation=90, va='center')
    plt.tight_layout(pad=1.2)

    plt.savefig(os.path.join(output_dir, "read_strand_types.png"), dpi=600, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, "read_strand_types.pdf"), bbox_inches='tight')

    handles, labels = plt.gca().get_legend_handles_labels()
    reordered_handles = []
    reordered_labels = []
    for i in range(3):
        reordered_handles.append(handles[i])
        reordered_handles.append(handles[i+3])
        reordered_labels.append(labels[i])
        reordered_labels.append(labels[i+3])
    fig_legend = plt.figure(figsize=(9, 1.2))
    fig_legend.legend(reordered_handles, reordered_labels, loc='center', ncol=3, fontsize=12, frameon=True, fancybox=True, handlelength=1.0)
    fig_legend.savefig(os.path.join(output_dir, "read_strand_types_legend.png"), dpi=600, bbox_inches='tight')
    fig_legend.savefig(os.path.join(output_dir, "read_strand_types_legend.pdf"), bbox_inches='tight')
    plt.close(fig_legend)
    plt.clf()

def plot_read_dist_types_gfeat2type(total_counts, output_dir='.'):
    plt.figure(figsize=(8, 8))
    plt.tight_layout()
    colors = ['#FBC2B5', '#FFA8A9', '#F786AA', '#A14A76',
                '#95F9E3', '#69EBD0', '#49D49D', '#558564']
    labels = ['Int Exonic', 'Int Intronic', 'Int Intergenic', 'Int Unassigned',
              'UMI Exonic', 'UMI Intronic', 'UMI Intergenic', 'UMI Unassigned']

    for i, row in enumerate(total_counts):
        sum_internal = sum(row[0:4])
        sum_umi = sum(row[4:8])
        row_internal = [x / sum_internal for x in row[0:4]]
        row_umi = [x / sum_umi for x in row[4:8]]
        total_counts[i] = [] + row_internal + row_umi

    data = np.array(total_counts)
    feature_names = ['Exon', 'Intron', 'Intergen', 'Unasgn']

    for i in range(8):
        avg = np.mean(data[:, i]) * 100
        print(f"{labels[i]}: {avg:.2f}%")

    fig, axs = plt.subplots(1, 4, figsize=(5, 5))
    axs = axs.flatten()
    for ax in axs:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    for feat_idx in range(4):
        ax = axs[feat_idx]
        feat_data = [data[:, feat_idx], data[:, feat_idx+4]]

        bp = ax.boxplot(feat_data, patch_artist=True, widths=0.7, showmeans=False,
                       medianprops=dict(linewidth=0))
        bp['boxes'][0].set_facecolor(colors[feat_idx])
        bp['boxes'][1].set_facecolor(colors[feat_idx+4])

        ax.set_xticks([1, 2])
        ax.set_xticklabels(['Internal', 'UMI'], fontsize=12)
        ax.tick_params(axis='x', pad=12)
        y_max = max([data[:, feat_idx].max(), data[:, feat_idx+4].max()])
        ceiling = np.ceil(y_max * 10 + 0.001) / 10

        if ceiling >= 0.6:
            tickes = [x * 0.2 for x in range(0, 6)]
            idx = min([i for i,x in enumerate(tickes) if x >= ceiling])
            ax.set_ylim(0, tickes[idx])
            ax.set_yticks(tickes[0: idx + 1])
            ax.set_yticklabels([f'{int(x*100)}%' for x in tickes][0: idx + 1], fontsize=12, rotation=90, va='center')
        elif ceiling >= 0.3:
            tickes = [x * 0.1 for x in range(0, 7)]
            idx = min([i for i,x in enumerate(tickes) if x >= ceiling])
            ax.set_ylim(0, tickes[idx])
            ax.set_yticks(tickes[0: idx + 1])
            ax.set_yticklabels([f'{int(x*100)}%' for x in tickes][0: idx + 1], fontsize=12, rotation=90, va='center')
        elif ceiling >= 0.2:
            tickes = [x * 0.05 for x in range(0, 7)]
            idx = min([i for i,x in enumerate(tickes) if x >= ceiling])
            ax.set_ylim(0, tickes[idx])
            ax.set_yticks(tickes[0: idx + 1])
            ax.set_yticklabels([f'{int(x*100)}%' for x in tickes][0: idx + 1], fontsize=12, rotation=90, va='center')
        else:
            tickes = [x * 0.05 for x in range(0, 5)]
            ax.set_ylim(0, 0.2)
            ax.set_yticks(tickes)
            ax.set_yticklabels([f'{int(x*100)}%' for x in tickes], fontsize=12, rotation=90, va='center')
        ax.set_xticks([])
        ax.set_title(feature_names[feat_idx], y=-0.075, fontsize=12)

    plt.subplots_adjust(wspace=0)
    plt.tight_layout(pad=2)
    plt.savefig(os.path.join(output_dir, "read_dist_feat2type_subplots.png"), dpi=600, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, "read_dist_feat2type_subplots.pdf"), bbox_inches='tight')

    from matplotlib.patches import Patch
    fig_legend = plt.figure(figsize=(2, 0.5))
    legend_elements = [Patch(facecolor='#95F9E3', label='UMI'),
                      Patch(facecolor='#FBC2B5', label='Internal')]
    fig_legend.legend(handles=legend_elements, loc='center', ncol=2, fontsize=12, frameon=True, fancybox=True, handlelength=1.0)
    fig_legend.savefig(os.path.join(output_dir, "read_dist_feat2type_subplots_legend.png"), dpi=600, bbox_inches='tight')
    fig_legend.savefig(os.path.join(output_dir, "read_dist_feat2type_subplots_legend.pdf"), bbox_inches='tight')
    plt.close(fig_legend)
    plt.clf()

def plot_strand_dist_strand2type(strand_counts, output_dir='.'):
    plt.figure(figsize=(8, 8))
    plt.tight_layout()
    colors = ['#FBC2B5', '#FFA8A9', '#A14A76',
                '#95F9E3', '#69EBD0', '#558564']
    labels = ['Internal FR', 'Internal RF', 'Internal Unknown',
              'UMI FR', 'UMI RF', 'UMI Unknown']

    for i, row in enumerate(strand_counts):
        sum_internal = sum(row[0:3])
        sum_umi = sum(row[3:6])
        assert sum_internal > 0.99 and sum_umi > 0.99
        assert sum_internal < 1.001 and sum_umi < 1.001
        row_internal = [x / sum_internal for x in row[0:3]]
        row_umi = [x / sum_umi for x in row[3:6]]
        strand_counts[i] = [] + row_internal + row_umi

    data = np.array(strand_counts)
    feature_names = ['FR', 'RF', 'Unkown']

    for i in range(6):
        avg = np.mean(data[:, i]) * 100
        print(f"{labels[i]}: {avg:.2f}%")

    fig, axs = plt.subplots(1, 3, figsize=(5, 5))
    axs = axs.flatten()
    for ax in axs:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    for feat_idx in range(3):
        ax = axs[feat_idx]
        feat_data = [data[:, feat_idx], data[:, feat_idx+3]]

        bp = ax.boxplot(feat_data, patch_artist=True, widths=0.7, showmeans=False,
                       medianprops=dict(linewidth=0))
        bp['boxes'][0].set_facecolor(colors[feat_idx])
        bp['boxes'][1].set_facecolor(colors[feat_idx+3])

        ax.set_xticks([1, 2])
        ax.set_xticklabels(['Internal', 'UMI'], fontsize=12)
        ax.tick_params(axis='x', pad=12)

        tickes = [x * 0.2 for x in range(0, 6)]
        ax.set_ylim(0, tickes[5])
        ax.set_yticks(tickes[0: 6])
        ax.set_yticklabels([f'{int(x*100)}%' for x in tickes][0: 6], fontsize=12, rotation=90, va='center')
        ax.set_xticks([])
        ax.set_title(feature_names[feat_idx], y=-0.075, fontsize=12)

    plt.subplots_adjust(wspace=0)
    plt.tight_layout(pad=2)
    plt.savefig(os.path.join(output_dir, "read_strand_strand2type_subplots.png"), dpi=600, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, "read_strand_strand2type_subplots.pdf"), bbox_inches='tight')

    from matplotlib.patches import Patch
    fig_legend = plt.figure(figsize=(2, 0.5))
    legend_elements = [Patch(facecolor='#95F9E3', label='UMI'),
                      Patch(facecolor='#FBC2B5', label='Internal')]
    fig_legend.legend(handles=legend_elements, loc='center', ncol=2, fontsize=12, frameon=True, fancybox=True, handlelength=1.0)
    fig_legend.savefig(os.path.join(output_dir, "read_strand_strand2type_subplots_legend.png"), dpi=600, bbox_inches='tight')
    fig_legend.savefig(os.path.join(output_dir, "read_strand_strand2type_subplots_legend.pdf"), bbox_inches='tight')
    plt.close(fig_legend)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--read-dist-dir', required=True,
                        help='Directory containing RSeQC read_distribution output files (*.dist)')
    parser.add_argument('--strand-dir', required=True,
                        help='Directory containing RSeQC infer_experiment output files (*.libtype.txt)')
    parser.add_argument('--output-dir', default='.',
                        help='Output directory for plots (default: current directory)')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Load read distribution files
    read_dist_files = []
    for root, _, files in os.walk(args.read_dist_dir):
        for file in files:
            if file.endswith('.dist'):
                read_dist_files.append(os.path.join(root, file))
    for f in read_dist_files:
        id = os.path.splitext(os.path.basename(f))[0]
        assert id.endswith(".internal") or id.endswith(".umied"), f"Unexpected file: {f}"
    read_dist_all_files = []
    for f in read_dist_files:
        id = os.path.splitext(os.path.basename(f))[0]
        read_dist = parse_read_dist(f)
        read_dist_all_files.append((id, read_dist))
    read_dist_all_files.sort(key=lambda x: x[0])

    # Load strand distribution files
    strand_files = []
    for root, _, files in os.walk(args.strand_dir):
        for file in files:
            if file.endswith('.libtype.txt'):
                strand_files.append(os.path.join(root, file))
    for f in strand_files:
        id = os.path.splitext(os.path.basename(f))[0]
        assert id.endswith(".internal.libtype") or id.endswith(".umied.libtype"), f"Unexpected file: {f}"
    strand_dist_all_files = []
    for f in strand_files:
        id = os.path.splitext(os.path.basename(f))[0]
        strand_dist = parse_libtype_file(f)
        strand_dist_all_files.append((id, strand_dist))
    strand_dist_all_files.sort(key=lambda x: x[0])

    assert len(read_dist_all_files) == len(strand_dist_all_files)

    read_counts = flatten_total_counts(read_dist_all_files, use_tags=False)
    plot_read_num_comp(deepcopy(read_counts), use_tags=False, output_dir=args.output_dir)

    read_counts = flatten_total_counts(read_dist_all_files, use_tags=True)
    plot_read_num_comp(deepcopy(read_counts), use_tags=True, output_dir=args.output_dir)

    read_dist = flatten_read_distributions(read_dist_all_files)
    plot_read_dist_types(deepcopy(read_dist), output_dir=args.output_dir)
    plot_read_dist_types_gfeat2type(deepcopy(read_dist), output_dir=args.output_dir)

    read_counts = flatten_total_counts(read_dist_all_files, use_tags=False)
    strand_counts = flatten_strand_counts(strand_dist_all_files)
    strand_dist = get_strand_dist(deepcopy(read_counts), deepcopy(strand_counts))
    plot_strand_dist(strand_dist, output_dir=args.output_dir)
    plot_strand_dist_strand2type(deepcopy(strand_counts), output_dir=args.output_dir)

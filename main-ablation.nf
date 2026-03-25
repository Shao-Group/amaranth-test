#!/usr/bin/env nextflow

// Ablation study workflow for Amaranth parameter tuning.
// Systematically varies ONE parameter at a time while keeping all others
// at permissive (no-effect) values to isolate each parameter's contribution.
//
// Usage: nextflow run main-ablation.nf -c params.config
// Set params.ablation_mode = "quick" for defaults-vs-permissive only (no sweeps).

def current_time = new Date().format('yyyy-MM-dd_HH-mm-ss')

// ============================================================================
// ABLATION STUDY
// Each experiment varies ONE parameter while keeping all others at permissive
// (no-effect) values, isolating each parameter's contribution.
// ============================================================================

process RUNAMARANTH_ABLATION {
    input:
    tuple val(id), path(bam), val(ablation_label), val(extra_args)

    output:
    tuple val(ablation_label), path("${id}.abl.${ablation_label}.gtf"), emit: gtf

    publishDir "${params.output_dir}/ablation/${ablation_label}-${current_time}", mode: 'copy', overwrite: true

    script:
    """
    amaranth -i ${bam} \
        ${extra_args} \
        -o ${id}.abl.${ablation_label}.gtf 
    """
}

process RUNGFFCOMPARE {
    input:
    tuple val(ablation_label), path(gtf)
    path reference_gtf

    output:
    tuple val(ablation_label), path("${gtf.baseName}.compare.stats"), emit: statsfile

    publishDir "${params.output_dir}/ablation-gffcpr/${ablation_label}-${current_time}", mode: 'copy', overwrite: true

    script:
    """
    gffcompare -M -N -r ${reference_gtf} -o ${gtf.baseName}.compare ${gtf}

    if [ -f ${gtf.baseName}.compare ]; then
        mv "${gtf.baseName}.compare" "${gtf.baseName}.compare.stats"
    fi
    """
}

process ABLATION_SUMMARY {
    input:
    val ablation_records

    output:
    path "ablation_summary.tsv", emit: summary

    publishDir "${params.output_dir}/ablation-summary-${current_time}", mode: 'copy', overwrite: true

    script:
    def header = "param_name\tparam_value\tamaranth_dir\tgffcompare_dir"
    def rows = ablation_records.collect { rec ->
        def label = rec
        def last_underscore = label.lastIndexOf('_')
        def pname = (last_underscore > 0) ? label.substring(0, last_underscore) : label
        def pval  = (last_underscore > 0) ? label.substring(last_underscore + 1) : label
        def amr_dir = "${params.output_dir}/ablation/${label}-${current_time}"
        def gff_dir = "${params.output_dir}/ablation-gffcpr/${label}-${current_time}"
        "${pname}\t${pval}\t${amr_dir}\t${gff_dir}"
    }.join('\n')
    """
    cat <<'ENDTSV' > ablation_summary.tsv
${header}
${rows}
ENDTSV
    """
}


workflow {
    println "Starting ablation workflow at ${current_time}"

    bams_ch = Channel.fromPath(params.bam_files)
        .map { bam ->
            def sample_id = bam.baseName
            tuple(sample_id, bam)
        }

    // --- Default values (tuned combination) ---
    def default_ir_part_v         = 0.5
    def default_ir_part_e         = 0.5
    def default_ir_full_v         = 1.0
    def default_ir_full_e         = 0.5
    def default_ir_full_i         = 10.0
    def default_umi_start_exon    = 1
    def default_umi_reads_bundle  = 1
    def default_umi_ratio_bundle  = 0.0

    // --- Permissive (no-effect) values ---
    // When sweeping one parameter, all others use these values
    def permissive_ir_part_v         = 0
    def permissive_ir_part_e         = 0
    def permissive_ir_full_v         = 0
    def permissive_ir_full_e         = 0
    def permissive_ir_full_i         = 999
    def permissive_umi_start_exon    = 0
    def permissive_umi_reads_bundle  = 0
    def permissive_umi_ratio_bundle  = 0.0

    // --- Ablation mode ---
    def ablation_mode = params.containsKey('ablation_mode') ? params.ablation_mode : "full"
    assert ablation_mode in ["full", "quick"] : "Invalid ablation_mode: ${ablation_mode}. Use 'full' or 'quick'."
    println "Ablation mode: ${ablation_mode}"

    // --- Sweep value lists ---
    def sweep_ir_part_v         = [0, 0.1, 0.5, 1.0, 2.0, 5.0]
    def sweep_ir_part_e         = [0, 0.1, 0.5, 1.0, 2.0, 5.0]
    def sweep_ir_full_v         = [0, 0.5, 1.0, 2.0, 5.0, 10.0]
    def sweep_ir_full_e         = [0, 0.1, 0.5, 1.0, 2.0, 5.0]
    def sweep_ir_full_i         = [0, 2.0, 5.0, 10.0, 20.0]
    def sweep_umi_start_exon    = [0, 1, 2, 5]
    def sweep_umi_reads_bundle  = [0, 1, 2, 5]
    def sweep_umi_ratio_bundle  = [0.0, 0.01, 0.05, 0.1]

    // --- Helper: build the extra_args string for one experiment ---
    def buildArgs = { ir_part_v, ir_part_e, ir_full_v, ir_full_e, ir_full_i,
                      umi_start_exon, umi_reads_bundle, umi_ratio_bundle ->
        return "--max-ir-part-ratio-v ${ir_part_v} " +
               "--max-ir-part-ratio-e ${ir_part_e} " +
               "--max-ir-full-ratio-v ${ir_full_v} " +
               "--max-ir-full-ratio-e ${ir_full_e} " +
               "--max-ir-full-ratio-i ${ir_full_i} " +
               "--min-umi-reads-start-exon ${umi_start_exon} " +
               "--min-umi-reads-bundle ${umi_reads_bundle} " +
               "--min-umi-ratio-bundle ${umi_ratio_bundle}"
    }

    // --- Generate ablation experiments ---
    def ablation_experiments = []

    // Baseline 1: all defaults (tuned combination)
    ablation_experiments << [
        "all_defaults",
        buildArgs(default_ir_part_v, default_ir_part_e, default_ir_full_v, default_ir_full_e,
                  default_ir_full_i, default_umi_start_exon,
                  default_umi_reads_bundle, default_umi_ratio_bundle)
    ]

    // Baseline 2: all permissive (no filtering)
    ablation_experiments << [
        "all_permissive",
        buildArgs(permissive_ir_part_v, permissive_ir_part_e, permissive_ir_full_v, permissive_ir_full_e,
                  permissive_ir_full_i, permissive_umi_start_exon,
                  permissive_umi_reads_bundle, permissive_umi_ratio_bundle)
    ]

    // --- Parameter sweeps (skipped in "quick" mode) ---
    if (ablation_mode == "full") {

    sweep_ir_part_v.each { v ->
        ablation_experiments << [
            "ir_part_v_${v}",
            buildArgs(v, permissive_ir_part_e, permissive_ir_full_v, permissive_ir_full_e,
                      permissive_ir_full_i, permissive_umi_start_exon,
                      permissive_umi_reads_bundle, permissive_umi_ratio_bundle)
        ]
    }

    sweep_ir_part_e.each { v ->
        ablation_experiments << [
            "ir_part_e_${v}",
            buildArgs(permissive_ir_part_v, v, permissive_ir_full_v, permissive_ir_full_e,
                      permissive_ir_full_i, permissive_umi_start_exon,
                      permissive_umi_reads_bundle, permissive_umi_ratio_bundle)
        ]
    }

    sweep_ir_full_v.each { v ->
        ablation_experiments << [
            "ir_full_v_${v}",
            buildArgs(permissive_ir_part_v, permissive_ir_part_e, v, permissive_ir_full_e,
                      permissive_ir_full_i, permissive_umi_start_exon,
                      permissive_umi_reads_bundle, permissive_umi_ratio_bundle)
        ]
    }

    sweep_ir_full_e.each { v ->
        ablation_experiments << [
            "ir_full_e_${v}",
            buildArgs(permissive_ir_part_v, permissive_ir_part_e, permissive_ir_full_v, v,
                      permissive_ir_full_i, permissive_umi_start_exon,
                      permissive_umi_reads_bundle, permissive_umi_ratio_bundle)
        ]
    }

    sweep_ir_full_i.each { v ->
        ablation_experiments << [
            "ir_full_i_${v}",
            buildArgs(permissive_ir_part_v, permissive_ir_part_e, permissive_ir_full_v, permissive_ir_full_e,
                      v, permissive_umi_start_exon,
                      permissive_umi_reads_bundle, permissive_umi_ratio_bundle)
        ]
    }

    sweep_umi_start_exon.each { v ->
        ablation_experiments << [
            "umi_start_exon_${v}",
            buildArgs(permissive_ir_part_v, permissive_ir_part_e, permissive_ir_full_v, permissive_ir_full_e,
                      permissive_ir_full_i, v,
                      permissive_umi_reads_bundle, permissive_umi_ratio_bundle)
        ]
    }

    sweep_umi_reads_bundle.each { v ->
        ablation_experiments << [
            "umi_reads_bundle_${v}",
            buildArgs(permissive_ir_part_v, permissive_ir_part_e, permissive_ir_full_v, permissive_ir_full_e,
                      permissive_ir_full_i, permissive_umi_start_exon,
                      v, permissive_umi_ratio_bundle) + " --both-umi-support true"
        ]
    }

    sweep_umi_ratio_bundle.each { v ->
        ablation_experiments << [
            "umi_ratio_bundle_${v}",
            buildArgs(permissive_ir_part_v, permissive_ir_part_e, permissive_ir_full_v, permissive_ir_full_e,
                      permissive_ir_full_i, permissive_umi_start_exon,
                      permissive_umi_reads_bundle, v) + " --both-umi-support true"
        ]
    }

    } // end if (ablation_mode == "full")

    println "Ablation: ${ablation_experiments.size()} experiments to run"

    ablation_ch = Channel.from(ablation_experiments)
        .map { item -> tuple(item[0], item[1]) }

    // Combine with bams: [id, bam, ablation_label, extra_args]
    ablation_input = bams_ch.combine(ablation_ch)

    // Run ablation
    abl_amr = RUNAMARANTH_ABLATION(ablation_input)

    // Ablation GFFCompare
    abl_gffcpr = RUNGFFCOMPARE(abl_amr.gtf, params.reference_gtf)

    // Collect all ablation labels for the summary table
    abl_labels = abl_gffcpr.statsfile
        .map { label, stats -> [label] }
        .collect()

    ABLATION_SUMMARY(abl_labels)
}

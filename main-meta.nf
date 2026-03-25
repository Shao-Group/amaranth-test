#!/usr/bin/env nextflow

// Meta-assembly transcript benchmarking workflow.
//
// Usage: nextflow run main-meta.nf -c params.config

def current_time = new Date().format('yyyy-MM-dd_HH-mm-ss')

process GETBEDFROMGTF {
    input:
    path gtf

    output:
    path "${gtf.baseName}.bed", emit: out

    publishDir "${params.output_dir}/ref_files", mode: 'copy', overwrite: true

    script:
    """
    gxf2bed -i ${gtf} -o "${gtf.baseName}.bed"
    """
}

process INFEREXPERIMENT {
    input:
    tuple val(id), path(bam), path(bed)

    output:
    tuple val(id), path(bam), path("${id}.libtype.txt"), emit: out

    publishDir "${params.output_dir}/libtype",
        mode: 'copy',
        overwrite: true,
        saveAs: { filename ->
            filename.endsWith(".libtype.txt") ? filename : null
        }

    script:
    """
    #!/bin/bash
    infer_experiment.py -r ${bed} -i ${bam} > ${id}.libtype.txt

    python -c '
    with open(\"${id}.libtype.txt\", \"r\") as f:
        for line in f:
            if \"failed to determine\" in line:
                pct_failed = float(line.split(\": \")[1])
            elif \"1++,1--,2+-,2-+\" in line:
                pct_FR = float(line.split(\": \")[1])
            elif \"1+-,1-+,2++,2--\" in line:
                pct_RF = float(line.split(\": \")[1])
        libtype = \"\"
        if pct_failed > 0.5:
            libtype = \"Unkown\"
        elif pct_FR >= pct_RF * 4 and pct_FR > 0.7:
            libtype = \"FR\"
        elif pct_RF >= pct_FR * 4 and pct_RF > 0.7:
            libtype = \"RF\"
        else:
            libtype = \"Unkown\"
        with open(\"${id}.libtype.txt\", \"a\") as f:
            f.write(libtype)
    '
    """
}

process RUNAMARANTH_META {
    input:
    path(bam_list)

    output:
    path "amaranth_meta_cell_*.gtf", emit: gtf

    publishDir "${params.output_dir}/amaranth-meta-${current_time}", mode: 'copy', overwrite: true

    script:
    """
    samtools merge -f merged.bam -b ${bam_list}
    
    mkdir meta
    amaranth --meta -i merged.bam \
        -o meta/amaranth 

    mkdir single
    for bam in \$(cat ${bam_list}); do
        amaranth -i \${bam} \
            -o single/\$(basename \${bam} .bam).amaranth
    done

    meta_files=(\$(ls meta/*.gtf | sort))
    single_files=(\$(ls single/*.gtf | sort))

    for i in \$(seq 0 \$((\${#meta_files[@]} - 1))); do
        gtfmerge union "\${meta_files[\$i]}" "\${single_files[\$i]}" > "amaranth_meta_cell_\${i}.gtf"
    done

    """
}

process RUNPSICLASS {
    input:
    path bam_list

    output:
    path "psiclass*.gtf", emit: gtf

    publishDir "${params.output_dir}/psiclass", mode: 'copy', overwrite: true

    script:
    def vd = params.containsKey('psiclass_vd') ? params.psiclass_vd : 1
    def threads = params.containsKey('psiclass_threads') ? params.psiclass_threads : 10
    """
    psiclass --lb ${bam_list} -o ./ --vd ${vd} -p ${threads}
    """
}

process RUNALETSCH {
    input:
    path bam_list

    output:
    path "*.gtf", emit: gtf

    publishDir "${params.output_dir}/aletsch", mode: 'copy', overwrite: true

    script:
    def c_param = params.containsKey('aletsch_c') ? params.aletsch_c : null
    """
    aletsch --profile -i ${bam_list} -p profile > profile.preview

    if [ -z "${c_param}" ]; then
        num=\$(wc -l < ${bam_list})
        c=\$((2 * num))
    else
        c=${c_param}
    fi

    aletsch -i ${bam_list} -o meta.gtf -d gtf -p profile -c \${c} > aletsch.log
    """
}

process RUNGFFCOMPARE {
    input:
    tuple path(gtf), val(toolname)
    path reference_gtf

    output:
    path "${gtf.baseName}.stats", emit: statsfile

    publishDir "${params.output_dir}/gffcpr-${current_time}", mode: 'copy', overwrite: true

    script:
    """
    gffcompare -M -N -r ${reference_gtf} -o ${gtf.baseName} ${gtf}

    if [ -f ${gtf.baseName} ]; then
        mv "${gtf.baseName}" "${gtf.baseName}.stats"
    fi

    # Assert stats file exists
    if [ ! -f "${gtf.baseName}.stats" ]; then
        echo "ERROR: GFFCompare output file "${gtf.baseName}.stats" not found"
        exit 1
    fi
    """
}

workflow {
    println "Starting meta-assembly workflow at ${current_time}"

    bams_ch = Channel.fromPath(params.bam_files)
        .map { bam ->
            def sample_id = bam.baseName
            tuple(sample_id, bam)
        }

    // Run tools
    bam_list = bams_ch.map { it[1] }.collectFile(name: 'bam.list', newLine: true)
    amr = RUNAMARANTH_META(bam_list)
    psi = RUNPSICLASS(bam_list)
    alt = RUNALETSCH(bam_list)

    // Channel gtf, run gffcompare
    gtf_assemblies_ch = amr.gtf.map { gtf -> tuple(gtf, 'amaranth-meta') }
        .mix(psi.gtf.map { gtf -> tuple(gtf, 'psiclass') })
        .mix(alt.gtf.map { gtf -> tuple(gtf, 'aletsch') })

    // evaluate gtf
    gffcpr = RUNGFFCOMPARE(gtf_assemblies_ch, params.reference_gtf)
}

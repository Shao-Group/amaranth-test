#!/usr/bin/env nextflow

// global variable
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
    def lib_type = ""
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
            libtype = \"Unkown\"     #libtype = \"Failed\"
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
    // Example:
    // This is PairEnd Data
    // Fraction of reads failed to determine: 0.0172
    // Fraction of reads explained by "1++,1--,2+-,2-+": 0.4903
    // Fraction of reads explained by "1+-,1-+,2++,2--": 0.4925
}

process RUNSTRINGTIE {
    input:
    tuple val(id), path(bam), path(libtype)

    output:
    path "${id}.stg.gtf", emit: gtf

    publishDir "${params.output_dir}/stringtie", mode: 'copy', overwrite: true

    script:
    """
    strand=\$(tail -n 1 ${libtype})
    echo "Strand: \${strand}"

    if [[ ! "\${strand}" =~ ^(RF|FR|Unkown)\$ ]]; then
        echo "Invalid strand: \${strand}"
        exit 1
    fi

    if [ "\${strand}" == "Unkown" ]; then
        lib_type=""
    elif [ "\${strand}" == "RF" ]; then
        lib_type="--rf"
    else
        lib_type="--fr"
    fi

    echo "stringtie version `stringtie --version`"
    stringtie -i ${bam} \${lib_type} -o ${id}.stg.gtf
    """
}

process RUNSCALLOP2 {
    input:
    tuple val(id), path(bam), path(libtype)

    output:
    path "${id}.sc2.gtf", emit: gtf

    publishDir "${params.output_dir}/scallop2", mode: 'copy', overwrite: true

    script:
    """
    strand=\$(tail -n 1 ${libtype})
    echo "Strand: \${strand}"

    if [[ ! "\${strand}" =~ ^(RF|FR|Unkown)\$ ]]; then
        echo "Invalid strand: \${strand}"
        exit 1
    fi

    if [ "\${strand}" == "Unkown" ]; then
        lib_type=""
    elif [ "\${strand}" == "RF" ]; then
        lib_type="--library_type first"
    else
        lib_type="--library_type second"
    fi
    
    echo "scallop2 version `scallop2 --version`"
    scallop2 -i ${bam} \${lib_type} -o ${id}.sc2.gtf
    """
}


process RUNAMARANTH {
    input:
    tuple val(id), path(bam), path(_libtype)

    output:
    path "${id}.amaranth.gtf", emit: gtf

    publishDir "${params.output_dir}/amaranth-${current_time}", mode: 'copy', overwrite: true

    script:

    """
    amaranth -i ${bam} \
                            --gene_name_prefix cell_${id}. \
                            -o ${id}.amaranth
    """
}

process RUNGFFCOMPARE {
    input:
    tuple path(gtf), val(toolname)
    path reference_gtf

    output:
    path "${gtf.baseName}.stats", emit: statsfile
    path "*.tmap", emit: tmapfile

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
    println "Starting workflow at ${current_time}"


    bams_ch = Channel.fromPath(params.bam_files)
        .map { bam -> 
            def sample_id = bam.baseName
            tuple(sample_id, bam)
        }

    bed = GETBEDFROMGTF(params.reference_gtf)

    lib = INFEREXPERIMENT(bams_ch.combine(bed.out))

    // Run tools
    sc2 = RUNSCALLOP2(lib)
    stg = RUNSTRINGTIE(lib)
    amr = RUNAMARANTH(lib)

    // Channel gtf, run IRtools
    gtf_assemblies_ch = sc2.gtf.map { gtf -> tuple(gtf, 'scallop2') }
        .mix(stg.gtf.map { gtf -> tuple(gtf, 'stringtie') })
        .mix(amr.gtf.map { gtf -> tuple(gtf, 'amaranth') })

    
    // evaluate gtf intron-chain level
    gffcpr = RUNGFFCOMPARE(gtf_assemblies_ch, params.reference_gtf)
}
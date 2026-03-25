// collect all stats files, compose tsv
process STATS2TSV {

    input:
    path statsfiles    // all stats files, from gffcpr.statsfile.collect()
    val  level         // intron chain or transcript
    val current_time  // current time for output directory

    output:
    path "*.tsv", emit: out

    publishDir "${params.output_dir}/tsv-${level}-${current_time}", mode: 'copy', overwrite: true

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    import sys
    import os
    from collections import defaultdict

    level_id_algo = defaultdict(lambda: defaultdict(lambda: dict()))

    for f in "${statsfiles}".split():
        with open(f, 'r') as file:
            fsplit = os.path.basename(f).split('.')
            fsplit[-1] == "stats"
            fsplit[-2] == "compare"
            id = fsplit[0]
            algo = '.'.join(fsplit[1:-2])
            lines = file.readlines()
            for line in lines:
                line = line.strip().lower()
                # precision
                if "level" in line:
                    line = line.replace('level', '')
                    level = line.split(':')[0].strip()          # e.g. 'Intron chain level'
                    entry = line.split(':')[1].strip().split()  # e.g. ['0.8', '|', '21.1', '|']
                    sensitivity = float(entry[0])
                    precision = float(entry[2])
                    level_id_algo[level + " sensitivity"][id][algo] = sensitivity
                    level_id_algo[level + " precision"][id][algo]   = precision
                # recall number
                if line.startswith("matching"):
                    line = line.replace('matching', '')
                    level = line.split(':')[0].strip()          # e.g. 'Matching intron chains'
                    level = level[:-1] if level.endswith('s') else level
                    recall = int(line.split(':')[1].strip())
                    level_id_algo[level + " number"][id][algo] = recall

    for level, id_algo in level_id_algo.items():
        df = pd.DataFrame(id_algo)
        df = df.sort_index()
        df['Average'] = df.mean(axis=1)
        df.insert(0, 'Average', df.pop('Average'))
        if level.endswith('number'):
            df = df.astype(int)
        else:
            df = df.round(2)  # precision 2 decimal
        df.to_csv(f"{'_'.join(level.split())}.tsv", sep='\\t')

    # Calculate ratio relative to 'sc2' row
    for level, id_algo in level_id_algo.items():
        df2 = pd.DataFrame(id_algo)
        df2 = df2.sort_index()
        for col in df2.columns:
            sc2_value = df2.loc['sc2', col]
            df2[col] = df2[col] / sc2_value
        df2['Average_FC'] = df2.mean(axis=1)
        df2.insert(0, 'Average_FC', df2.pop('Average_FC'))
        df2 = df2.round(2)
        df2.to_csv(f"fc_{'_'.join(level.split())}.tsv", sep='\\t')
    """
}


// format tsv, if qsv is installed
process FORMATTSV {
    input:
    path tsv
    val  level         // intron chain or transcript
    val current_time  // current time for output directory

    output:
    path "q.*.tsv", emit: out_q, optional: true
    path "t.*.tsv", emit: out_t, optional: true

    publishDir "${params.output_dir}/qsv-${level}-${current_time}", mode: 'copy', overwrite: true

    script:
    """
    if ! command -v qsv &> /dev/null; then
        echo "qsv is required but not installed; Skipping formatting"
    else
        qsv table ${tsv} > q.${tsv.baseName}.tsv
        qsv transpose ${tsv} > t.${tsv.baseName}.tsv
    fi
    """
}

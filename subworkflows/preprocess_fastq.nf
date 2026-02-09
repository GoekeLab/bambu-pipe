process QFILTER{
    container "ghcr.io/ch99l/bambu-pipe:latest"
    label "low_cpu"
    label "low_mem"
    label "long"

    input:
    tuple val(sample), path(fastq), val(meta)

    output:
    tuple val(sample), path("${sample}_qfilter.fastq*"), val(meta)

    script:
    """
    pigz -p $task.cpus -d -c $fastq | chopper -q $params.qfilter_threshold -t $task.cpus > ${sample}_qfilter.fastq

    if [[ $params.compress_intermediate == "true" ]]; then
        pigz -p $task.cpus ${sample}_qfilter.fastq
    fi
    """
}

process DEMULTIPLEX{  
	container "ghcr.io/ch99l/bambu-pipe:latest"
    label "medium_cpu"
    label "low_mem"
    label "long"

	input: 
	tuple val(sample), path(fastq), val(meta)

	output:
    tuple val(sample), path("${sample}_flexiplexfilter_reads.fastq*"), val(meta)

	script:
	"""
    if [[ $params.compress_intermediate == "true" ]]; then
        pigz -p $task.cpus -d -c $fastq > flexiplex_input_reads.fastq
        fastq_file="flexiplex_input_reads.fastq"
    else
        fastq_file="$fastq"
    fi

    whitelist=\$(awk -F',' -v chem=$meta.chemistry '\$1 == chem {print \$2}' ${projectDir}/10x_config/barcode_config.csv)
    if [[ \$whitelist == *.gz ]]; then
        pigz -p $task.cpus -d -c \$whitelist > whitelist.txt
    else
        cp \$whitelist whitelist.txt
    fi

    IFS=',' read -r _ left_flank barcode umi right_flank \
    < <(awk -F',' -v chem=$meta.chemistry '\$1 == chem' ${projectDir}/10x_config/flank_seq_config.csv)
    flank_seq="-x \$left_flank -b \$barcode -u \$umi -x \$right_flank"

	flexiplex -p $task.cpus \$flank_seq -f 0 \$fastq_file
	flexiplex-filter -w whitelist.txt --outfile my_filtered_barcode_list.txt flexiplex_barcodes_counts.txt 
    flexiplex -p $task.cpus -k my_filtered_barcode_list.txt \$flank_seq -f $params.flexiplex_f -e $params.flexiplex_e \$fastq_file > ${sample}_flexiplexfilter_reads.fastq

    if [[ $params.compress_intermediate == "true" ]]; then
        rm flexiplex_input_reads.fastq
        pigz -p $task.cpus ${sample}_flexiplexfilter_reads.fastq
    fi
    """
}

process TRIM_AND_ORIENT{
    container "ghcr.io/ch99l/bambu-pipe:latest"
    label "medium_cpu"
    label "low_mem"
    label "short"

    input:
    tuple val(sample), path(fastq), val(meta)
    
    output:
    tuple val(sample), path("${sample}_preprocessed_reads.fastq*"), val(meta)

    script:
    """
    IFS=',' read -r _ fwd_primer_f fwd_primer_r rev_primer_f rev_primer_r TSO_f TSO_r \
    < <(awk -F',' -v chem=$meta.chemistry '\$1 == chem' ${projectDir}/10x_config/adapter_seq_config.csv)

    if [[ $meta.chemistry == 10x5* ]]; then
        cutadapt -a \$rev_primer_f --cores $task.cpus $fastq | \
        cutadapt -b \$fwd_primer_f -b \$fwd_primer_r -b \$TSO_f -b \$TSO_r -b \$rev_primer_f -b \$rev_primer_r --action none --discard \
        --cores $task.cpus -o ${sample}_preprocessed_reads.fastq - # For 5' preparation kits, reads are already in the transcript direction 

    elif [[ $meta.chemistry == 10x3* || $meta.chemistry == visium-v* ]]; then
        cutadapt -a \$rev_primer_f --cores $task.cpus $fastq | \
        cutadapt -b \$fwd_primer_f -b \$fwd_primer_r -b \$rev_primer_f -b \$rev_primer_r --action none --discard --cores $task.cpus - | \
        reverse_complement_fastq.py -i - -o ${sample}_preprocessed_reads.fastq # For 3' preparation kits, orient reads in the transcript direction to improve minimap alignment
        
    fi

    if [[ $params.compress_intermediate == "true" ]]; then
        pigz -p $task.cpus ${sample}_preprocessed_reads.fastq
    fi
    """
}

workflow PREPROCESS_FASTQ {
    take: 
    ch_fastq

    main:
    // process fastq samples
    ch_flexiplex_in = params.qscore_filtering ? QFILTER(ch_fastq) : ch_fastq // skip quality score filtering if params.qscore_filtering is false
    ch_flexiplex_in | DEMULTIPLEX | TRIM_AND_ORIENT // chain preprocessing steps

    emit:
    TRIM_AND_ORIENT.out
}
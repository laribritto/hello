/*
 * pipeline input parameters
 */
params.reads = "$projectDir/data/ggal/gut_{1,2}.fq"
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
params.multiqc = "$projectDir/multiqc"
params.outdir = "resultados"

// Validação dos parâmetros obrigatórios
if (params.reads == null || params.transcriptome_file == null || params.outdir == null) {
    log.error "Parâmetros obrigatórios não foram fornecidos corretamente. Verifique se os caminhos para os arquivos de leitura (reads), o arquivo de transcriptoma (transcriptome_file) e o diretório de saída (outdir) estão definidos."
    System.exit(1) // Encerra o script com código de erro 1
}

log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    transcriptome: ${params.transcriptome_file}
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    """
    .stripIndent()

/*
 * define the `index` process that creates a binary index
 * given the transcriptome file
 */
process INDEX {
    container 'quay.io/biocontainers/salmon:1.10.1--hecfa306_2'
    conda 'bioconda::salmon=1.10.1'
    input:
    path transcriptome

    output:
    path 'salmon_index'

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i salmon_index
    """
}

process QUANTIFICATION {
    container 'quay.io/biocontainers/salmon:1.10.1--hecfa306_2'
    conda 'bioconda::salmon=1.10.1'
    // Processo para realizar a quantificação com Salmon
    tag "Salmon on $sample_id"
    publishDir params.outdir, mode:'copy'

    input:
    path salmon_index
    tuple val(sample_id), path(reads)

    output:
    path "$sample_id"

    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $salmon_index -1 ${reads[0]} -2 ${reads[1]} -o $sample_id
    """
}

process FASTQC {
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    conda 'bioconda::fastqc=0.11.8'

    tag "FASTQC on $sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}

process MULTIQC {
    container 'quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0'
    conda 'bioconda::fastqc=1.14'
    publishDir params.outdir, mode:'copy'

    input:
    path '*'

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}

workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    index_ch = INDEX(params.transcriptome_file)
    quant_ch = QUANTIFICATION(index_ch, read_pairs_ch)
    fastqc_ch = FASTQC(read_pairs_ch)
    MULTIQC(quant_ch.mix(fastqc_ch).collect())
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}

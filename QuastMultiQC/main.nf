#!/usr/bin/env nextflow

// Created by Gerrit Botha (grbot), edited by me (BournSupremacy)
// Get sample info from sample sheet
// Minimum information that is needed in the sample sheet are SampleID and fasta file location of that sample
Channel.fromPath( file(params.sample_sheet) )
        .splitCsv(header: true, sep: '\t')
        .map{row ->
            def sample_id = row['SampleID']
            def fasta_file = file(row['fasta'])
            return [ sample_id, fasta_file ]
        }.set{samples}

out_dir = file(params.out_dir)
ref = file(params.ref)
nt = params.nt

process getQuast{
    tag { "${params.project_name}.${sample_id}.gQ" }
    cpus { "${nt}" }
    publishDir "${out_dir}/${sample_id}", mode: 'copy', overwrite: false
   
    input:
	  set val(sample_id), file(fasta_file) from samples

    output:
	  set val(sample_id), file ("*") into quast_files

    script:
    """
    quast -r ${ref} \
    -o ${sample_id}.quast \
    --fragmented --eukaryote \
    --conserved-genes-finding \
    -t ${task.cpus} \
    ${fasta_file}
    """
}

process runMultiQC{
    tag { "${params.project_name}.rMQC" }
    publishDir "${out_dir}/", mode: 'copy', overwrite: false

    input:
        file('*') from quast_files.collect()
	
    output:
        file('multiqc_report.html')

    """
    multiqc .
    """
}

workflow.onComplete {

    println ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )
}

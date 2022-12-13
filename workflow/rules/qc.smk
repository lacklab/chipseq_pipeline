rule fastqc:
    input:
        "links/{raw}_{run}.fastq.gz"
    output:
        "qc/{raw}.{run}_fastqc.zip"
    shell:
        """
        fastqc {input} -o qc
        """

rule multiqc:
    input:
        get_multiqc
    output:
        "qc/multiqc_report.html"
    shell:
        """
        cd qc/ && multiqc .
        """
# Rule: Prefetch SRA data
#rule sra_prefetch:
#    output:
#        temp("sra-data/{SRA}/{SRA}.sra")
#    conda:
#        "../envs/sra.yaml"
#    shell:
#        """
#        prefetch -O sra-data {wildcards.SRA}
#        """

# Rule: Convert SRA to FASTQ using parallel-fastq-dump
rule parallel_fastq_dump:
    output:
        r1 = "sra-data/{srr}_1.fastq.gz",
        r2 = "sra-data/{srr}_2.fastq.gz"
    conda:
        "../envs/sra.yaml"
    threads:
        64
    shell:
        r"""
        LIBTYPE=$(awk -v srr={wildcards.srr} '$5==srr{{print $4}}' config/samples.tsv)

        if [ "$LIBTYPE" = "Single" ]; then
            parallel-fastq-dump -t {threads} --split-files --gzip -s {wildcards.srr} -O sra-data
            touch {output.r2}  # create dummy r2 for Single-end
        elif [ "$LIBTYPE" = "Paired" ]; then
            parallel-fastq-dump -t {threads} --split-files --gzip -s {wildcards.srr} -O sra-data
        else
            echo "Error: Unknown library type for {wildcards.srr}" >&2
            exit 1
        fi
        """
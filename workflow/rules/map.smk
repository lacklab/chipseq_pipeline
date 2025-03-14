# Rule: Map reads to the reference genome using BWA

rule map_bwa:
    input:
        trimmed_fq1="trimmed/{raw}_1.trimmed.fastq.gz",
        trimmed_fq2="trimmed/{raw}_2.trimmed.fastq.gz"
    output:
        raw=temp("results_{ref}/mapping/{raw}.raw.bam"),  # Temporary raw BAM file
        log="qc/bwa/{ref}:{raw}.bwa.log"
    params:
        idx=lambda wildcards: references[wildcards.ref]["BWA_IDX"],
        lib=lambda wildcards: get_lib(wildcards)
    threads:
        16  # Number of threads to use
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        if [[ "{params.lib}" == "Single" ]]; then
            bwa mem -t {threads} {params.idx} {input.trimmed_fq1} 2> {output.log} \
            | samtools view -bS - > {output.raw}
        elif [[ "{params.lib}" == "Paired" ]]; then
            bwa mem -t {threads} {params.idx} {input} 2> {output.log} \
            | samtools view -bS - > {output.raw}
        fi
        """

# Rule: Process BAM file (coordinate sorting, fixing mates, marking duplicates)
rule bam_process:
    input:
        "results_{ref}/mapping/{raw}.raw.bam"  # Raw BAM file from mapping
    output:
        temp("results_{ref}/mapping/{raw}.coorsorted.bam")  # Coordinate-sorted BAM file
    params:
        config["OUTPUT"]["BAMPROCESS_PARAMS"]  # Parameters for BAM processing
    threads:
        16
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        samtools view -h {params} {input} \
        | samtools fixmate -m -@ {threads} - - \
        | samtools sort -@ {threads} -m 10G - \
        | samtools markdup -@ {threads} - {output}
        """

# Rule: Filter BAM file (remove mitochondrial reads, apply blacklist, etc.)
rule bam_filter:
    input:
        "results_{ref}/mapping/{raw}.coorsorted.bam"  # Coordinate-sorted BAM file
    output:
        temp("results_{ref}/mapping/{raw}.filtered.bam")  # Filtered BAM file
    params:
        fa = lambda wildcards: references[wildcards.ref]["FA"],  # Reference genome FASTA
        p  = lambda wildcards: "-F 3852 -f 2" if get_lib(wildcards) == "Paired" else "-F 3852",  # Filtering parameters (e.g., flags)
        bl = lambda wildcards: references[wildcards.ref]["BLACKLIST"]  # Path to blacklist file
    threads:
        16
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        samtools view {input} | egrep -v "chrM" \
        | samtools view -b -@ {threads} -T {params.fa} {params.p} \
        | bedtools intersect -nonamecheck -v -abam stdin -b {params.bl} > {output}
        """

# Rule: Merge replicates into a final BAM file
rule bam_merge:
    input:
        get_units
    output:
        "results_{ref}/mapping/{name}.final.bam"  # Final merged BAM file
    threads:
        16
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        if [[ $(echo {input} | grep ' ') ]]; then
            samtools merge -@ {threads} -o {output} {input}
        else
            mv {input} {output}
        fi
        samtools index {output}
        """

# Rule: Create pseudoreplicates from the final BAM file
# Will be deprecated!
rule pseudoreps:
    input:
        "results_{ref}/mapping/{name}.final.bam"  # Final BAM file
    output:
        pr1 = "results_{ref}/mapping/{name}.pr1.bam",  # First pseudoreplicate
        pr2 = "results_{ref}/mapping/{name}.pr2.bam"   # Second pseudoreplicate
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        samtools index {input}
        samtools view -b --subsample 0.5 {input} > {output.pr1}
        samtools view -b --subsample 0.5 {input} > {output.pr2}
        """
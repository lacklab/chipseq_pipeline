rule genomecov:
    input:
        bam = "results_{ref}/mapping/{name}.final.bam"
    output:
        bg = temp("results_{ref}/bigwig/{name}.genomecov.{norm}.bg"),
        bw = "results_{ref}/bigwig/{name}.genomecov.{norm}.bw",
        params_file = temp("results_{ref}/bigwig/{name}.params.{norm}.txt")  # Stores BAM parameters
    params:
        chrSizes = lambda wildcards: references[wildcards.ref]["CHROM_SIZES"]
    threads:
        4
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        # Run the Python script to determine BAM processing parameters
        python workflow/scripts/determine_bam_params.py --bam {input.bam} --norm {wildcards.norm} --output {output.params_file}

        # Generate genome coverage tracks
        bedtools genomecov -bg -ibam {input.bam} $(cat {output.params_file}) | \
        sort -k1,1 -k2,2n --parallel={threads} > {output.bg}
        
        bedGraphToBigWig {output.bg} {params.chrSizes} {output.bw}
        """

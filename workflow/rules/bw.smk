# Rule: Generate genome coverage tracks in BedGraph and BigWig formats
rule genomecov:
    input:
        bam = "results_{ref}/mapping/{name}.final.bam"  # Input final BAM file
    output:
        bg = temp("results_{ref}/bigwig/{name}.genomecov.{norm}.bg"),  # Temporary BedGraph file
        bw = "results_{ref}/bigwig/{name}.genomecov.{norm}.bw"  # Output BigWig file
    params:
        chrSizes = lambda wildcards: references[wildcards.ref]["CHROM_SIZES"]  # Path to chromosome sizes file
    threads:
        4  # Number of threads to use
    run:
        # Open the BAM file to determine read type (paired or single-end) and calculate scaling
        bam = pysam.AlignmentFile(input["bam"])
        first_read = next(bam)
        if first_read.is_paired:
            # Paired-end read parameters
            pc = " -pc "  # Use proper pairs
            fl = " "  # No fragment size adjustment for paired-end
            scale = (f" -scale {10**6 / sum([(r.flag & 64) == 64 for r in bam])}"
                     if wildcards.norm == "FPM" else " ")  # Scale by FPM if required
        else:
            # Single-end read parameters
            pc = " "  # No proper pairs for single-end
            fl = f" -fs {len(first_read.seq)}"  # Set fragment size to read length
            scale = f" -scale {10**6 / bam.mapped}" if wildcards.norm == "FPM" else " "  # Scale by FPM if required
        bam.close()

        # Generate genome coverage tracks
        shell("""
            bedtools genomecov \
            -bg -ibam {input.bam} {pc} {fl} {scale} | \
            sort -k1,1 -k2,2n --parallel={threads} > {output.bg}
            
            bedGraphToBigWig {output.bg} {params.chrSizes} {output.bw}
        """)
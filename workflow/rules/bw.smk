rule genomecov:
    input:
        bam = "results_{ref}/mapping/{raw}.final.bam",
    output:
        bg = temp("results_{ref}/bigwig/{raw}.genomecov.{norm}.bg"),
        bw = "results_{ref}/bigwig/{raw}.genomecov.{norm}.bw"
    params:
        chrSizes = config["REF"]["CHROM_SIZES"]
    threads:
        16
    run:
        bam = pysam.AlignmentFile(input["bam"])
        a = next(bam)
        fl = " " if a.is_paired else f" -fs {len(a.seq)}"
        pc = " -pc " if a.is_paired else " "
        scale = f" -scale {10**6 / sum([(a.flag & 64) == 64 for a in bam])}" if wildcards.norm == "RPM" else " "
        shell("""
            bedtools genomecov \
            -bg -ibam {input.bam} {pc} {fl} {scale} \
            | sort -k1,1 -k2,2n --parallel={threads} > {output.bg}
            bedGraphToBigWig {output.bg} {params.chrSizes} {output.bw}
            """)
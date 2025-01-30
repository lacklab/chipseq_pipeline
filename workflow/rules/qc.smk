# Rule: Generate alignment statistics using samtools stats
rule stats:
    input:
        "results_{ref}/mapping/{raw}.bam"  # Input BAM file
    output:
        "qc/samtools/stats/{ref}:{raw}.stats"  # Output stats file
    params:
        fa = lambda wildcards: references[wildcards.ref]["FA"]  # Reference genome path
    conda:
        "../envs/qc.yaml"
    shell:
        """
        samtools stats --reference {params.fa} {input} > {output}
        """

# Rule: Generate flagstat statistics using samtools flagstat
rule flagstat:
    input:
        "results_{ref}/mapping/{raw}.bam"  # Input BAM file
    output:
        "qc/samtools/flagstat/{ref}:{raw}.flagstat"  # Output flagstat file
    conda:
        "../envs/qc.yaml"
    shell:
        """
        samtools flagstat {input} > {output}
        """

# Rule: Generate idxstats statistics using samtools idxstats
rule idxstats:
    input:
        "results_{ref}/mapping/{raw}.bam"  # Input BAM file
    output:
        "qc/samtools/idxstats/{ref}:{raw}.idxstats"  # Output idxstats file
    conda:
        "../envs/qc.yaml"
    shell:
        """
        samtools idxstats {input} > {output}
        """

# Rule: Move MACS2 peak calling QC output
rule macs_qc:
    input:
        "results_{ref}/peaks/{name}_{q}_peaks.xls"  # Input MACS2 peaks file
    output:
        "qc/macs/{ref}:{name}_{q}_peaks.xls"  # QC directory for MACS2 peaks file
    conda:
        "../envs/qc.yaml"
    shell:
        """
        ln -s $(readlink -f {input}) $(readlink -f {output})
        """

# Rule: Annotate peaks and generate a summary for QC
from collections import Counter
rule annotatepeaks_qc:
    input:
        "results_{ref}/annot/{name}_{q}_annotatepeaks.txt"  # Input annotated peaks file
    output:
        "qc/homer/{ref}:{name}_{q}_summary_mqc.txt"  # Summary QC file
    run:
        header = ["INTERGENIC", "INTRON ", "PROMOTER-TSS ", "EXON ", "3' UTR ", "5' UTR ", "TTS ", "NON-CODING "]
        with open(output[0], "w") as f:
            f.write(assets["annotatepeaks"])  # Pre-written asset for annotatepeaks
            tmp = pd.read_table(input[0])  # Read annotation file
            if tmp.shape[0] == 0:
                nAnnot = dict(zip(header, [0] * len(header)))  # No annotations found
            else:
                tmp["shortAnn"] = tmp["Annotation"].str.split("(", expand=True)[0].str.upper()
                nAnnot = Counter(tmp["shortAnn"])  # Count annotations
            for k in header:
                f.write(f"{k}\t{nAnnot.get(k, 0)}\n")  # Write counts to summary file

# Rule: Compute FRiP (Fraction of Reads in Peaks) values


rule frip:
    input:
        bams=expand("results_{{ref}}/mapping/{name}.final.bam", name=samples.loc[samples['Control'] != '-', "Name"].tolist()),   # Fetch BAM files for FRIP calculation
        peak=expand("results_{{ref}}/peaks/{name}_{q}_peaks.narrowPeak", name=samples.loc[samples['Control'] != '-', "Name"].tolist(), q=config['OUTPUT']['MACS_THRESHOLD'])    # Fetch peak files for FRIP calculation
    output:
        "qc/{ref}:frip_mqc.tsv"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        python workflow/scripts/calculate_frip.py \
        --bams {input.bams} \
        --peaks {input.peak} \
        --output {output}
        """

# Rule: Generate a MultiQC report
rule multiqc:
    input:
        get_multiqc  # Collect all files for MultiQC
    output:
        "qc/multiqc_report.html"  # MultiQC output report
    conda:
        "../envs/qc.yaml"
    shell:
        """
        cd qc/ && multiqc .
        """

# TODO:
# - Integrate FRiP values into 'multiqc_data/multiqc_general_stats.txt'.
# - Add custom content to the MultiQC report (https://multiqc.info/docs/#custom-content).
# - Consider output examples like those from nf-core's ChIP-seq pipeline.
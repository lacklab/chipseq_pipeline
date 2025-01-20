# Rule: Annotate peaks using HOMER's annotatePeaks.pl
rule homer_annotatepeaks:
    input:
        # Input narrowPeak file
        "results_{ref}/peaks/{name}_{q}_peaks.narrowPeak"
    output:
        # Output annotated peaks file
        "results_{ref}/annot/{name}_{q}_annotatepeaks.txt"
    shell:
        """
        annotatePeaks.pl {input} {wildcards.ref} > {output}
        """
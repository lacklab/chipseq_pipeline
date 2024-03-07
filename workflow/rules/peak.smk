rule macs:
    input:
        get_macs_i
    output:
        "results_{ref}/peaks/{raw}_{p}_{q}_peaks.xls",
        "results_{ref}/peaks/{raw}_{p}_{q}_peaks.narrowPeak"
    threads:
        16
    params:
        get_macs_p
    shell:
        """
        macs3 callpeak \
            {params}  \
            -g hs -n results_{wildcards.ref}/peaks/{wildcards.raw}_{wildcards.p}_{wildcards.q} -B -q {wildcards.q}
        """

#config['OUTPUT']['IDR_THRESHOLD']
rule idr:
    input:
        unpack(get_idr_i)
    output:
        "results_{ref}/idr/{raw}_{q}_idr.narrowPeak"
    params:
        q = config['OUTPUT']['IDR_THRESHOLD']
    shell:
        """
        ~/miniconda/envs/notebook/bin/idr \
            --samples {input.pr1} {input.pr2} \
            --peak-list {input.final} \
            --input-file-type narrowPeak --output-file {output} \
            --rank signal.value --soft-idr-threshold {params.q} \
            --plot --use-best-multisummit-IDR
        """

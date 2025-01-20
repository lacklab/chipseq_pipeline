# Rule: Call peaks using MACS3
rule macs:
    input:
        get_macs_i  # Function to retrieve MACS input files
    output:
        # MACS3 output files: Excel summary and narrowPeak file
        peaks_xls="results_{ref}/peaks/{name}_{q}_peaks.xls",
        narrow_peak="results_{ref}/peaks/{name}_{q}_peaks.narrowPeak"
    threads:
        4  # Number of threads to use
    params:
        get_macs_p  # Function to retrieve MACS parameters
    shell:
        """
        macs3 callpeak \
            {params} \
            -g hs \
            -n results_{wildcards.ref}/peaks/{wildcards.name}_{wildcards.q} \
            -q {wildcards.q}  # FDR threshold for peak calling
        """

# Rule: Perform IDR (Irreproducible Discovery Rate) analysis
# DEPARCATED
#rule idr:
#    input:
#        unpack(get_idr_i)  # Unpack pseudoreplicates and final peak files
#    output:
#        idr_peak="results_{ref}/idr/{raw}_{q}_idr.narrowPeak"  # IDR-filtered peaks file
#    params:
#        q = config['OUTPUT']['IDR_THRESHOLD']  # IDR threshold from config
#    shell:
#        """
#        ~/miniconda/envs/notebook/bin/idr \
#            --samples {input.pr1} {input.pr2} \  # Pseudoreplicate files
#            --peak-list {input.final} \  # Final peaks file
#            --input-file-type narrowPeak \  # Input file format
#            --output-file {output.idr_peak} \  # Output file
#            --rank signal.value \  # Rank peaks by signal value
#            --soft-idr-threshold {params.q} \  # IDR threshold
#            --plot \  # Generate diagnostic plots
#            --use-best-multisummit-IDR  # Use best summit for multisummit peaks
#        """
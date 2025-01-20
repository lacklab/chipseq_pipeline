# Rule: Download BED files from ChIP-Atlas
rule ChipAtlasBed:
    output:
        # Output BED file for the specified SRX and threshold
        "results_{ref}/cabeds/srx/{srx}.{threshold}.bed"
    threads:
        16  # Number of threads to use for downloading
    shell:
        """
        url=https://chip-atlas.dbcls.jp/data/{wildcards.ref}/eachData/bed{wildcards.threshold}/{wildcards.srx}.{wildcards.threshold}.bed
        if wget $url -O- &>/dev/null; then
            /home/ualtintas/apps/aria2-1.35.0/src/aria2c -x {threads} -s {threads} \
                $url \
                -o {output}
        else
            touch {output}  # Create an empty file if the download fails
        fi
        """

# Rule: Link ChIP-Atlas BED files to a common directory
rule ChipAtlasBed_link:
    input:
        lambda wildcards: f"results_{wildcards.ref}/cabeds/srx/{assets['gsm2srx'][wildcards.gsm]}.{wildcards.threshold}.bed"  # Function to retrieve BED file paths
    output:
        # Output symlink to the specified BED file
        "results_{ref}/cabeds/{raw}.{gsm}.{threshold}.bed"
    shell:
        """
        ln -s $(readlink -f {input}) {output}
        """

# Rule: Download BigWig files from ChIP-Atlas
rule ChipAtlasBigwig:
    output:
        # Output BigWig file for the specified SRX
        "results_{ref}/cabws/srx/{srx}.bw"
    threads:
        16  # Number of threads to use for downloading
    shell:
        """
        url=https://chip-atlas.dbcls.jp/data/{wildcards.ref}/eachData/bw/{wildcards.srx}.bw
        if wget $url -O- &>/dev/null; then
            /home/ualtintas/apps/aria2-1.35.0/src/aria2c -x {threads} -s {threads} \
                $url \
                -o {output}
        else
            touch {output}  # Create an empty file if the download fails
        fi
        """

# Rule: Link ChIP-Atlas BigWig files to a common directory
rule ChipAtlasBigwig_link:
    input:
        lambda wildcards: f"results_{wildcards.ref}/cabws/srx/{assets['gsm2srx'][wildcards.gsm]}.bw" # Function to retrieve BigWig file paths
    output:
        # Output symlink to the specified BigWig file
        "results_{ref}/cabws/{raw}.{gsm}.bw"
    shell:
        """
        ln -s $(readlink -f {input}) {output}
        """
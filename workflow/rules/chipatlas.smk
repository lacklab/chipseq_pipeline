rule ChipAtlasBed:
    output:
        "results_{ref}/cabeds/srx/{srx}.{threshold}.bed"
    threads:
        16
    shell:
        """
        url=https://chip-atlas.dbcls.jp/data/{wildcards.ref}/eachData/bed{wildcards.threshold}/{wildcards.srx}.{wildcards.threshold}.bed
        if [[ $(wget $url -O-) ]] 2>/dev/null;
            then
                /home/ualtintas/apps/aria2-1.35.0/src/aria2c -x {threads} -s {threads} \
                $url \
                -o {output}
            else
                touch {output}
        fi
        """

rule ChipAtlasBed_link:
    input:
        get_cabeds
    output:
        "results_{ref}/cabeds/{raw}.{gsm}.{threshold}.bed"
    shell:
        """
        ln -s `readlink -f {input}` {output}
        """



rule ChipAtlasBigwig:
    output:
        "results_{ref}/cabws/srx/{srx}.bw"
    threads:
        16
    shell:
        """
        url=https://chip-atlas.dbcls.jp/data/{wildcards.ref}/eachData/bw/{wildcards.srx}.bw
        if [[ $(wget $url -O-) ]] 2>/dev/null;
            then
                /home/ualtintas/apps/aria2-1.35.0/src/aria2c -x {threads} -s {threads} \
                $url \
                -o {output}
            else
                touch {output}
        fi
        """

rule ChipAtlasBigwig_link:
    input:
        get_cabws
    output:
        "results_{ref}/cabws/{raw}.{gsm}.bw"
    shell:
        """
        ln -s `readlink -f {input}` {output}
        """
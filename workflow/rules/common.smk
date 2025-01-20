import os
import urllib.request
import pandas as pd
import pysam

# Load configuration and samples
# Uncomment and use if `yaml` and `Path` are available
import yaml
from pathlib import Path
references = yaml.safe_load(Path("config/references.yaml").read_text())

# `samples` includes 'IP/CONTROL' pairs
samples = pd.read_table(config["SAMPLES"])
samples["Raw"] = samples["Name"] + "_" + samples["Unit"].astype(str)

# Asset management
assets = {}
assets_path = "assets"

# Load annotatepeaks asset
annotatepeaks_file = os.path.join(assets_path, "annotatepeaks.asset")
with open(annotatepeaks_file, "r") as f:
    assets["annotatepeaks"] = "".join(f.readlines())

# Download or update experiment list
experiment_list_url = "https://chip-atlas.dbcls.jp/data/metadata/experimentList.tab"
experiment_list_file = os.path.join(assets_path, "experimentList.tab")
if not os.path.exists(experiment_list_file) or config["UPDATE_CHIPATLAS"]:
    urllib.request.urlretrieve(experiment_list_url, experiment_list_file)

# Parse experiment list for GSM and SRX mappings
gsm2srx, srx2gsm = [], []
with open(experiment_list_file, "r") as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 9 and "SRX" in parts[0] and parts[8] != '-':
            gsm = parts[8].split(':')[0]
            srx = parts[0]
            gsm2srx.append((gsm, srx))
            srx2gsm.append((srx, gsm))

assets["gsm2srx"] = dict(gsm2srx)
assets["srx2gsm"] = dict(srx2gsm)

# Utility functions
def get_lib(wildcards):
    return samples.loc[samples["Raw"] == wildcards.raw, "Library"].unique()[0]

def get_fastq(wildcards, fastq_col):
    return samples.loc[samples["Raw"] == wildcards.raw, fastq_col].unique()[0]


def get_group_deeptools(wildcards):
    return expand("results_{{ref}}/mapping/{name}.final.bam", name=samples["Name"].tolist())

def get_units(wildcards):
    return expand("results_{{ref}}/mapping/{rep}.filtered.bam", rep=samples.loc[samples["Name"] == wildcards.name, "Raw"].unique()) 



##########################

def get_multiqc(wildcards):
    out = []
    
    # List of common quality control file types and tools
    qc_tools_1 = {
        "fastqc": [
            "{raw}_1_fastqc.html", 
            "{raw}_2_fastqc.html"
        ],
        "trimgalore": [
            "{raw}_1.fastq.gz_trimming_report.txt",
            "{raw}_2.fastq.gz_trimming_report.txt",
            "{raw}_1.trimmed_fastqc.html",
            "{raw}_2.trimmed_fastqc.html"
        ],
        "samtools": [
            "flagstat/{ref}:{raw}.coorsorted.flagstat",
            "idxstats/{ref}:{raw}.coorsorted.idxstats",
            "stats/{ref}:{raw}.coorsorted.stats"
        ],
        "bwa": [
            "{ref}:{raw}.bwa.log"
        ],
        "deeptools": [
            "{ref}:all_bam.bamSummary.npz",
            "{ref}:all_bam.plotCorrelation.mat.tab",
            "{ref}:all_bam.plotFingerprint.qcmetrics.txt",
            "{ref}:all_bam.plotFingerprint.raw.txt",
            "{ref}:all_bam.plotPCA.tab"
        ]
    }
    qc_tools_2 = {
         "samtools": [
            "flagstat/{ref}:{name}.final.flagstat",
            "idxstats/{ref}:{name}.final.idxstats",
            "stats/{ref}:{name}.final.stats"
        ],
        "macs": [
            "{ref}:{name}_{q}_peaks.xls"
        ],
        "homer":
        [
            "{ref}:{name}_{q}_summary_mqc.txt"
        ] # TODO: duplications
    }
    # Iterate through each sample and append all files based on the defined templates
    for ref in config['OUTPUT']['REF']:
        for q in config['OUTPUT']['MACS_THRESHOLD']:
            for _, row in samples.iterrows():
                raw = row['Raw']
                name = row['Name']

                is_control = row['Control'] == '-'
                
                # Generate output paths for each tool and file pattern
                for tool, patterns in qc_tools_1.items():
                    for pattern in patterns:
                        out.append(f"qc/{tool}/{pattern.format(raw=raw, ref=ref, q=q)}")

                for tool, patterns in qc_tools_2.items():
                    for pattern in patterns:
                        if is_control:
                            continue
                        out.append(f"qc/{tool}/{pattern.format(name=name, ref=ref, q=q)}")
   
    # Add FRIP score file (outside the loop as a single file)
        out.append(f"qc/{ref}:frip_mqc.tsv")

    
    return expand(out)


# map.smk functions
def get_fqs(wildcards):
    fq1 = get_fastq(wildcards, "Fastq1")
    lib = get_lib(wildcards)
    if "SRR" in fq1:
        return (f"sra-data/{fq1}_1.fastq.gz", f"sra-data/{fq1}_2.fastq.gz") if lib == "Paired" else f"sra-data/{fq1}_1.fastq.gz"
    fq2 = get_fastq(wildcards, "Fastq2")
    return (fq1, fq2) if lib == "Paired" else fq1

# peak.smk functions
def get_macs_i(wildcards):
    control = samples.loc[samples["Name"] == wildcards.name, "Control"].unique()[0]
    inputs = [f"results_{wildcards.ref}/mapping/{wildcards.name}.final.bam"]
    if control != "-":
        inputs.append(f"results_{wildcards.ref}/mapping/{control}.final.bam")
    return inputs


def get_macs_p(wildcards):
    lib = samples.loc[samples["Name"] == wildcards.name, "Library"].unique()[0]
    input_c = samples.loc[samples["Name"] == wildcards.name, "Control"].unique()[0]
    params = f"-t results_{wildcards.ref}/mapping/{wildcards.name}.final.bam "
    if input_c != "-":
        params += f"-c results_{wildcards.ref}/mapping/{input_c}.final.bam "
    params += "-f BAMPE" if lib == "Paired" else "-f BAM"
    return params

# Outputs
bwNorm = config["OUTPUT"]["BW_NORMALIZATIONS"]
outputs = []

if config["OUTPUT"]["RUN"]["QC"]:
    outputs.append("qc/multiqc_report.html")

if config["OUTPUT"]["RUN"]["PEAKS"]:
    outputs += [
        f"results_{ref}/peaks/{row['Name']}_{q}_peaks.narrowPeak"
        for ref in config['OUTPUT']['REF']
        for q in config['OUTPUT']['MACS_THRESHOLD']
        for i, row in samples.iterrows()
		if row['Control'] != '-'

    ]

if config["OUTPUT"]["RUN"]["BWS"]:
    outputs += [
        f"results_{ref}/bigwig/{raw}.genomecov.{norm}.bw"
        for ref in config['OUTPUT']['REF']
        for raw in samples["Name"]
        for norm in bwNorm
    ]

if config["OUTPUT"]["RUN"]["CHIPATLASBED"]:
    caq = config["OUTPUT"]["CHIPATLASBED_THRESHOLD"]
    outputs += [
        f"results_{ref}/cabeds/{samples.loc[samples['GSM'] == gsm, 'Name'].iloc[0]}.{gsm}.{caq}.bed"
        for ref in config['OUTPUT']['REF']
        for gsm in samples["GSM"].unique()
        if gsm in assets["gsm2srx"]
    ]

if config["OUTPUT"]["RUN"]["CHIPATLASBIGWIG"]:
    outputs += [
        f"results_{ref}/cabws/{samples.loc[samples['GSM'] == gsm, 'Name'].iloc[0]}.{gsm}.bw"
        for ref in config['OUTPUT']['REF']
        for gsm in samples["GSM"].unique()
        if gsm in assets["gsm2srx"]
    ]


print(outputs)
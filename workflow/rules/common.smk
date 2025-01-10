import pandas as pd
import urllib.request
import os.path
import pysam
#import yaml
#from pathlib import Path
#config = yaml.safe_load(Path("config/config.yaml").read_text())

# `samples` includes 'IP/CONTROL' pairs 


samples = pd.read_table(config["SAMPLES"])
samples["Raw"] = samples["Name"] + "_" + samples["Unit"].astype(str)

#samples = samples.loc[:26,:]


#samples = pd.read_table(config["SAMPLES"])
#if config['OUTPUT']['RUN']['CHIPATLASBED'] or config['OUTPUT']['RUN']['CHIPATLASBIGWIG']:
#	geo_samples = pd.read_table(config["GEOSAMPLES"])

# `units` includes all files that need to be preprocessed
#units = pd.read_table(config["UNITS"])
#units["Raw"] = units["Name"] + "_" + units["Unit"].astype(str)


# >>> Assets >>>
assets = {}
with open("assets/annotatepeaks.asset", "r") as f:
	assets["annotatepeaks"]	= ""
	for line in f.readlines():
		assets["annotatepeaks"]	+= line
		
if not os.path.exists("assets/experimentList.tab"):
	urllib.request.urlretrieve("https://chip-atlas.dbcls.jp/data/metadata/experimentList.tab", "assets/experimentList.tab")
elif config["UPDATE_CHIPATLAS"]:
	urllib.request.urlretrieve("https://chip-atlas.dbcls.jp/data/metadata/experimentList.tab", "assets/experimentList.tab")

gsm2srx = []
srx2gsm = []
with open("assets/experimentList.tab", "r") as f:
    for line in f.readlines():
        line = line.split('\t')
        if len(line) < 9:
            print(line)
        if (line[0].find("SRX") == -1) or (line[8] == '-'):
            continue
        
        gsm2srx.append((line[8].split(':')[0], line[0]))
        srx2gsm.append((line[0], line[8].split(':')[0]))

assets["srx2gsm"] = dict(srx2gsm)
assets["gsm2srx"] = dict(gsm2srx)
# <<< Assets <<<


ref = config["OUTPUT"]["REF"]
q = config["OUTPUT"]["MACS_THRESHOLD"]
# >>> utils >>>
def get_lib(wildcards):
	return samples.loc[samples["Raw"] == wildcards.raw, "Library"].unique()[0]

def get_units(wildcards):
	return samples.loc[samples["Name"] == wildcards.raw, "Raw"].unique()

def get_fq1(wildcards):
	return samples.loc[samples["Raw"] == wildcards.raw, "Fastq1"].unique()[0]

def get_fq2(wildcards):
	return samples.loc[samples["Raw"] == wildcards.raw, "Fastq2"].unique()[0]

def get_contol(wildcards):
	return samples.loc[samples["Name"] == wildcards.raw, "Control"].unique()[0]
# <<< utils <<<

# >>> `qc.smk` >>>
fastqc_map = {}
for u in samples["Fastq1"].tolist() + samples["Fastq2"].tolist():
	if (u.find("gz") != -1) or (u.find("zip") != -1):
		fastqc_map[u.rsplit(".", 2)[0].split("/")[-1]] = u
	elif u.find("SRR") != -1:
		fastqc_map[u+"_1"] = f"sra-data/{u}_1.fastq.gz"
		fastqc_map[u+"_2"] = f"sra-data/{u}_2.fastq.gz"
	elif u != "-":
		fastqc_map[u.rsplit(".", 1)[0].split("/")[-1]] = u

def get_fastqc(wildcards):
	return fastqc_map[wildcards.raw]

def get_annotatepeaks(wildcards):
	out = []
	for i, row in samples.iterrows():
		out = list(set(out.append(f"qc/{ref}:{row['Name']}.annotatePeaks.txt")))
	return expand(out)

def get_frip_b(wildcards):
	out = []
	for i, row in samples.iterrows():
		out = list(set(out.append(f"results_{ref}/mapping/{row['Name']}.final.bam")))
	return out

def get_frip_p(wildcards):
	out = []
	for i, row in samples.iterrows():
		out = list(set(out.append(f"results_{ref}/peaks/{row['Name']}_{q}_peaks.narrowPeak")))
	return expand(out)

def get_multiqc(wildcards):
	out = []
	for i, row in samples.iterrows():
		lib = row["Library"]
		fq1 = row["Fastq1"].split("/")[-1]
		fq2 = row["Fastq2"].split("/")[-1]
		if fq1.find("SRR") != -1:
			ext1 = "_1"
			ext2 = "_2"
			fq2 = fq1
		else:
			ext1 = ""
			ext2 = ""
		if fq1.find("gz") != -1:
			fq1 = fq1.rsplit(".", 2)[0]
			fq2 = fq2.rsplit(".", 2)[0]
		if lib == "Single":
			out.append(f"qc/fastqc/{fq1+ext1}_fastqc.zip")
		elif lib == "Paired":
			out.append(f"qc/fastqc/{fq1+ext1}_fastqc.zip")
			out.append(f"qc/fastqc/{fq2+ext2}_fastqc.zip")
		out.append(f"qc/flagstats/{ref}:{row['Raw']}.raw")
		out.append(f"qc/flagstats/{ref}:{row['Name']}.final")
		out.append(f"qc/stats/{ref}:{row['Name']}.final")
	for i, row in samples.iterrows():
		out.append(f"qc/annot/{ref}:{row['Name']}_{q}.summary_mqc.txt")
		out.append(f"qc/macs/{ref}:{row['Name']}_{q}_peaks.xls")
	out.append("qc/frip_mqc.tsv")
	return expand(out)
# <<< `qc.smk` <<<

# >>> `map.smk` functions >>>
def get_fqs(wildcards):
	#check_out = checkpoints.parallel_fastq_dump.get(**wildcards).output
#
	#srr = wildcards.srr
#
#
	#name = samples['Name']
	#unit = samples['Unit']



	#name, unit = wildcards.raw.rsplit("_",1)
	fq1 = samples.loc[samples["Raw"] == wildcards.raw,"Fastq1"].unique()[0]
	source = str(fq1).find("SRR") != -1
	print(source, 'hey')
	lib = get_lib(wildcards)
	if source:
		srr = fq1
		if lib == "Single":
			return f"sra-data/{srr}_1.fastq.gz"
		elif lib == "Paired":
			return f"sra-data/{srr}_1.fastq.gz", f"sra-data/{srr}_2.fastq.gz"
	else:
		fq1 = get_fq1(wildcards)
		fq2 = get_fq2(wildcards)
		if lib == "Single":
			return fq1
		elif lib == "Paired":
			return fq1, fq2

def get_filter_p(wildcards):
	lib = get_lib(wildcards)
	if lib == "Single":
		return "-F 3852"
	elif lib == "Paired":
		return "-F 3852 -f 2"

def get_reps(wildcards):
	reps = get_units(wildcards)
	return expand("results_{{ref}}/mapping/{rep}.filtered.bam", rep=reps)
# <<< `map.smk` functions <<<

# >>> `peak.smk` functions >>>
def get_macs_p(wildcards):
	lib = samples.loc[samples["Name"] == wildcards.raw, "Library"].unique()[0]
	input_c = get_contol(wildcards)
	param = f"-t results_{wildcards.ref}/mapping/{wildcards.raw}.{wildcards.p}.bam "
	if not input_c == "-":
		param += f" -c results_{wildcards.ref}/mapping/{input_c}.final.bam "
	if lib == "Single":
		param += " -f BAM"
	elif lib == "Paired":
		param += " -f BAMPE"
	return param

def get_macs_i(wildcards):
	inputs = []
	input_c = get_contol(wildcards)
	if not input_c == "-":
		inputs.append(f"results_{wildcards.ref}/mapping/{input_c}.final.bam")
	inputs.append(f"results_{wildcards.ref}/mapping/{wildcards.raw}.{wildcards.p}.bam")
	return inputs

def get_idr_i(wildcards):
		qm = config['OUTPUT']['MACS_THRESHOLD']
		#qi = config['OUTPUT']['IDR_THRESHOLD']
		return {
			'pr1' : f"results_{wildcards.ref}/peaks/{wildcards.raw}_pr1_{qm}_peaks.narrowPeak",
        	'pr2' : f"results_{wildcards.ref}/peaks/{wildcards.raw}_pr2_{qm}_peaks.narrowPeak",
        	'final' : f"results_{wildcards.ref}/peaks/{wildcards.raw}_final_{qm}_peaks.narrowPeak"
			}




# <<< `peak.smk` functions <<<


# >>> `chipatlas.smk` functions >>>
def get_cabeds(wildcards):
	#gsm = units.loc[units["Name"] == wildcards.raw, "GSM"].unique()[0]
	srx = assets["gsm2srx"][wildcards.gsm]
	return f"results_{wildcards.ref}/cabeds/srx/{srx}.{wildcards.threshold}.bed"


def get_cabws(wildcards):
	#gsm = units.loc[units["Name"] == wildcards.raw, "GSM"].unique()[0]
	srx = assets["gsm2srx"][wildcards.gsm]
	return f"results_{wildcards.ref}/cabws/srx/{srx}.bw"
# <<< `chipatlas.smk` functions <<<

# >>> OUTPUTS >>>

bwNorm = config["OUTPUT"]["BW_NORMALIZATIONS"]

outputs = []

if config["OUTPUT"]["RUN"]["QC"]:
	outputs += ["qc/multiqc_report.html"]

qi = config['OUTPUT']['IDR_THRESHOLD']
if config["OUTPUT"]["RUN"]["PEAKS"]:
	outputs += [
		f"results_{ref}/idr/{raw}_{qi}_idr.narrowPeak"
		for raw in samples["Name"]
		if raw.find('input') == -1
	]

if config["OUTPUT"]["RUN"]["BWS"]:
	outputs += [
		f"results_{ref}/bigwig/{raw}.genomecov.{norm}.bw"
		for raw in samples["Name"]
		for norm in bwNorm
	]

gsm2n = dict(zip(samples['GSM'], samples['Name']))
caq = config['OUTPUT']['CHIPATLASBED_THRESHOLD']
if config["OUTPUT"]["RUN"]["CHIPATLASBED"]:
	outputs += [
		f"results_{ref}/cabeds/{gsm2n[gsm]}.{gsm}.{caq}.bed"
		for gsm in samples["GSM"].unique()
		if gsm in assets["gsm2srx"]
	]

samples["GSM"].unique()

if config["OUTPUT"]["RUN"]["CHIPATLASBIGWIG"]:
	outputs += [
		f"results_{ref}/cabws/{gsm2n[gsm]}.{gsm}.bw"
		for gsm in samples["GSM"].unique()
		if gsm in assets["gsm2srx"]
	]

# <<< OUTPUTS <<<
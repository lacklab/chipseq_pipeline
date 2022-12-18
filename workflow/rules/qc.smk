rule fastqc:
	input:
		get_fastqc
	output:
		"qc/{raw}_fastqc.zip"
	shell:
		"""
		fastqc {input} -o qc
		"""

rule stats:
	input:
		"results_{ref}/mapping/{raw}.bam"
	output:
		"qc/{ref}:{raw}.bam_stats"
	shell:
		"""
		samtools stats {input} > {output}
		"""

rule flagstat:
	input:
		"results_{ref}/mapping/{raw}.bam"
	output:
		"qc/{ref}:{raw}.bam_flagstat"
	shell:
		"""
		samtools flagstat {input} > {output}
		"""

rule homer_annotatepeaks:
	input:
		"results_{ref}/peaks/{raw}_{q}_peaks.narrowPeak"
	output:
		"qc/{ref}:{raw}_{q}.annotatepeaks.txt"
	shell:
		"""
		annotatePeaks.pl {input} {ref} > {output}
		"""

from collections import Counter
import re
rule annotatepeaks_merge:
	input:
		get_annotatepeaks
	output:
		"qc/annotatepeaks.summary_mqc.tsv"
	run:
		header = ["intron","intergenic", "promoter-tss", "exon", "tts", "5' utr", "3' utr"]
		with open(output) as f:
			f.write("\t".join(["Sample Name"] + header))
			for s in input:
				sample_name = re.split(s, "/|\.")[1]
				tmp = pd.read_table(s)
				tmp = tmp.rename(columns={tmp.columns[0]: "PeakID"})
				tmp["shortAnn"] = tmp["Annotation"].str.split("(", expand=True)[0].str.upper()
				nAnnot = Counter(tmp["shortAnn"])
				f.write("\t".join([sample_name] + [nAnnot[h] for h in header]))
				
# TODO: This can be written in a script.	

rule multiqc:
	input:
		get_multiqc
	output:
		"qc/multiqc_report.html"
	shell:	
		"""
		cd qc/ && multiqc .
		"""

# TODO: get FRiP values and add to 'multiqc_data/multiqc_general_stats.txt'
# TODO: work on custom-content (https://multiqc.info/docs/#custom-content)
# TODO: QC output can include outputs; such as https://nf-co.re/chipseq/1.2.2/output
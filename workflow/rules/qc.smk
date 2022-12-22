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

from collections import Counter
rule annotatepeaks_qc:
	input:
		"results_{ref}/annot/{raw}_{q}.annotatepeaks.txt"
	output:
		"qc/{ref}:{raw}_{q}.annotatepeaks.summary_mqc.txt"
	run:
		with open(output[0], "w") as f:
			f.write(assets["annotatepeaks"])
			tmp = pd.read_table(input).rename(columns={tmp.columns[0]: "PeakID"})
			tmp["shortAnn"] = tmp["Annotation"].str.split("(", expand=True)[0].str.upper()
			nAnnot = Counter(tmp["shortAnn"])
			for k, v in nAnnot.items():
				f.write(f"{k}\t{v}\n")

import deeptools.countReadsPerBin as crpb
import pysam		
rule frip:
	input:
		get_frip
	output:
		"qc/frip_mqc.tsv"
	run:
		with open(output[0], "w") as f:
			f.write("# plot_type: 'generalstats'\n")
			f.write("Sample Name\tFRiP")
			for b, p in zip(input["BAMS"], input["PEAKS"]):
				cr = crpb.CountReadsPerBin([b],bedFile=[p],numberOfProcessors=10)
				reads_at_peaks = cr.run()
				total = reads_at_peaks.sum(axis=0)
				bam = pysam.AlignmentFile(b)
				name = b.split("/")[-1].split("_peaks")[0]
				f.write(f"{name}\t{float(total[0]) / bam.mapped}")

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
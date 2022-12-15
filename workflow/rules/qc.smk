rule fastqc:
	input:
		get_fastqc
	output:
		"qc/{raw}_fastqc.zip"
	shell:
		"""
		fastqc {input} -o qc
		"""

rule multiqc:
	input:
		get_multiqc
	output:
		"qc/multiqc_report.html"
	shell:	
		"""
		cd qc/ && multiqc .
		"""        
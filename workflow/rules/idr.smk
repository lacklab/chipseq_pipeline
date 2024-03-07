

# TF Chip-seq 2: https://www.encodeproject.org/pipelines/ENCPL367MAS/ 


# ===================
# Create tagAlign file
# ===================
module add bedtools/2.26.0

# Create SE tagAlign file
FINAL_TA_FILE="${FINAL_BAM_PREFIX}.SE.tagAlign.gz"

bedtools bamtobed -i ${FINAL_BAM_FILE} \
| awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' \
| gzip -nc > ${FINAL_TA_FILE}


#########

# TAGALIGN PAIRED

# ===================
# Create tagAlign file
# ===================
module add bedtools/2.26.0


# ================
# Create BEDPE file
# ================
FINAL_BEDPE_FILE="${FINAL_NMSRT_BAM_PREFIX}.bedpe.gz"
bedtools bamtobed -bedpe -mate1 -i ${FINAL_NMSRT_BAM_FILE} | gzip -nc > ${FINAL_BEDPE_FILE}

zcat ${FINAL_BEDPE_FILE} \
| awk 'BEGIN{OFS="\t"}{\
	printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10\
	}' \
| gzip -nc > ${FINAL_TA_FILE}




#############

# SINGLE

PR_PREFIX="${OFPREFIX}.filt.nodup"
PR1_TA_FILE="${PR_PREFIX}.SE.pr1.tagAlign.gz"
PR2_TA_FILE="${PR_PREFIX}.SE.pr2.tagAlign.gz"

# Get total number of read pairs
nlines=$( zcat ${FINAL_TA_FILE} | wc -l )
nlines=$(( (nlines + 1) / 2 ))

# Shuffle and split BED file into 2 equal parts
zcat ${FINAL_TA_FILE} | shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:$(zcat -f ${FINAL_TA_FILE} | wc -c) -nosalt </dev/zero 2>/dev/null)  | split -d -l ${nlines} - ${PR_PREFIX} # Will produce ${PR_PREFIX}00 and ${PR_PREFIX}01

# Convert reads into standard tagAlign file
gzip -nc “${PR_PREFIX}00" > ${PR1_TA_FILE}
rm "${PR_PREFIX}00"
gzip -nc “${PR_PREFIX}01" > ${PR2_TA_FILE}
rm "${PR_PREFIX}01"




###########

# PAIRED

# ========================
# Create pseudoReplicates
# =======================
PR_PREFIX="${OFPREFIX}.filt.nodup"
PR1_TA_FILE="${PR_PREFIX}.PE2SE.pr1.tagAlign.gz"
PR2_TA_FILE="${PR_PREFIX}.PE2SE.pr2.tagAlign.gz"
joined=”temp.bedpe”

# Make temporary fake BEDPE file from FINAL_TA_FILE
zcat ${FINAL_TA_FILE} | sed 'N;s/\n/\t/' > $joined

# Get total number of read pairs
nlines=$( zcat ${joined} | wc -l )
nlines=$(( (nlines + 1) / 2 ))

# Shuffle and split BEDPE file into 2 equal parts
zcat ${joined} | shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:$(zcat -f ${FINAL_TA_FILE} | wc -c) -nosalt </dev/zero 2>/dev/null)  | split -d -l ${nlines} - ${PR_PREFIX} # Will produce ${PR_PREFIX}00 and ${PR_PREFIX}01

# Convert fake BEDPE into standard tagAlign file
awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12} "${PR_PREFIX}00" | gzip -nc > ${PR1_TA_FILE}
rm "${PR_PREFIX}00"
awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}'"${PR_PREFIX}01" | gzip -nc > ${PR2_TA_FILE}
rm "${PR_PREFIX}01"
rm -f ${joined}


############

# POOLED

# ========================
# Create pooled datasets
# =======================
REP1_TA_FILE=”${DATASET_PREFIX}.Rep1.tagAlign.gz”
REP2_TA_FILE=”${DATASET_PREFIX}.Rep2.tagAlign.gz”
POOLED_TA_FILE=”${DATASET_PREFIX}.Rep0.tagAlign.gz”

zcat ${REP1_TA_FILE} ${REP2_TA_FILE} | gzip -nc > ${POOLED_TA_FILE}

# ========================
# Create pooled pseudoreplicates
# =======================
REP1_PR1_TA_FILE=”${DATASET_PREFIX}.Rep1.pr1.tagAlign.gz”
REP1_PR2_TA_FILE=”${DATASET_PREFIX}.Rep1.pr2.tagAlign.gz”

REP2_PR1_TA_FILE=”${DATASET_PREFIX}.Rep2.pr1.tagAlign.gz”
REP2_PR2_TA_FILE=”${DATASET_PREFIX}.Rep2.pr2.tagAlign.gz”

PPR1_TA_FILE=”${DATASET_PREFIX}.Rep0.pr1.tagAlign.gz”
PPR2_TA_FILE=”${DATASET_PREFIX}.Rep0.pr2.tagAlign.gz”

zcat ${REP1_PR1_TA_FILE} ${REP2_PR1_TA_FILE} | gzip -nc > ${PPR1_TA_FILE}
zcat ${REP1_PR2_TA_FILE} ${REP2_PR2_TA_FILE} | gzip -nc > ${PPR2_TA_FILE}






# encode3 pipeline: https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit

# SPP usage: https://github.com/ENCODE-DCC/chip-seq-pipeline2/issues/109



IDR_THRESH=0.05

# =============================
# Perform IDR analysis.
# Generate a plot and IDR output with additional columns including IDR scores.
# =============================

idr \
--samples ${REP1_PEAK_FILE} ${REP2_PEAK_FILE} \
--peak-list ${POOLED_PEAK_FILE} \
--input-file-type narrowPeak --output-file ${IDR_OUTPUT} \
--rank signal.value --soft-idr-threshold ${IDR_THRESH} \
--plot --use-best-multisummit-IDR


# =============================
# Get peaks passing IDR threshold of 5%
# =============================

# Col 12: -log10(global IDR value)

IDR_THRESH_TRANSFORMED=$(awk -v p=${IDR_THRESH} 'BEGIN{print -log(p)/log(10)}')

awk \
'BEGIN{OFS="\t"} \
$12>='"${IDR_THRESH_TRANSFORMED}"' \
{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' \
${IDR_OUTPUT} \
| sort | uniq \
| sort -k7n,7n \
| gzip -nc > ${REP1_VS_REP2}.IDR0.05.narrowPeak.gz

NPEAKS_IDR=$(zcat ${REP1_VS_REP2}.IDR0.05.narrowPeak.gz | wc -l)


# =============================
# Filter using black list
# =============================
bedtools intersect -v -a ${REP1_VS_REP2}.IDR0.05.narrowPeak.gz -b ${BLACKLIST} | grep -P 'chr[\dXY]+[ \t]' | awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}’
| gzip -nc > ${REP1_VS_REP2}.IDR0.05.filt.narrowPeak.gz

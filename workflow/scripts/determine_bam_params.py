import argparse
import pysam
import os

def determine_bam_params(bam_file, norm):
    """Determine whether the BAM file is paired-end or single-end, and compute necessary parameters."""
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        first_read = next(bam)
        is_paired = first_read.is_paired

        if is_paired:
            pc = "-pc"
            fl = ""
            total_reads = sum(1 for r in bam if (r.flag & 64) == 64)  # Count first mates only
        else:
            pc = ""
            fl = "-fs " + str(len(first_read.seq))
            total_reads = bam.mapped  # Use total mapped reads

        # Compute scaling factor for FPM normalization
        scale = f"-scale {10**6 / total_reads}" if norm == "FPM" else ""

    return pc, fl, scale

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Determine BAM processing parameters")
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--norm", required=True, help="Normalization method (FPM or none)")
    parser.add_argument("--output", required=True, help="Output file to store parameters")

    args = parser.parse_args()
    
    pc, fl, scale = determine_bam_params(args.bam, args.norm)

    # Write parameters to a file (to be read by Snakemake)
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, "w") as f:
        f.write(f"{pc} {fl} {scale}")

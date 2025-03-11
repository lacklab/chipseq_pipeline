import argparse
import deeptools.countReadsPerBin as crpb
import pysam
import numpy as np

def calculate_frip(bam_files, peak_files, output_file):
    with open(output_file, "w") as f:
        f.write("# plot_type: 'generalstats'\n")
        f.write("Sample Name\tFRiP\tNumber of Peaks\tMedian Fragment Length\n")

        for bam_path, peak_path in zip(bam_files, peak_files):
            # Count reads at peaks using deepTools
            cr = crpb.CountReadsPerBin([bam_path], bedFile=[peak_path], numberOfProcessors=10)
            reads_at_peaks = cr.run()
            total_reads_at_peaks = reads_at_peaks.sum(axis=0)

            # Calculate total mapped reads using pysam
            bam = pysam.AlignmentFile(bam_path)
            total_mapped_reads = bam.mapped

            # Count number of peaks
            with open(peak_path, 'r') as peak_file:
                num_peaks = sum(1 for _ in peak_file)

            # Calculate median fragment length
            fragment_lengths = [
                abs(read.template_length) for read in bam.fetch() if read.is_proper_pair
            ]
            median_fragment_length = np.median(fragment_lengths)

            # Compute FRIP score
            sample_name = peak_path.split("/")[-1].split("_peaks")[0]
            frip_score = float(total_reads_at_peaks[0]) / total_mapped_reads

            # Write results
            f.write(f"{sample_name}\t{frip_score:.4f}\t{num_peaks}\t{median_fragment_length:.2f}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate FRIP scores for given BAM and peak files.")
    parser.add_argument("--bams", nargs="+", required=True, help="List of BAM files")
    parser.add_argument("--peaks", nargs="+", required=True, help="List of peak files")
    parser.add_argument("--output", required=True, help="Output TSV file")
    args = parser.parse_args()

    calculate_frip(args.bams, args.peaks, args.output)

import argparse
import os
import matplotlib.pyplot as plt
import subprocess
import pandas as pd
import pysam

def parse_args():
    """define args from the command line"""
    parser = argparse.ArgumentParser(description="Map reads to a reference genome and make diagnostic plots to assess assembly quality.")
    parser.add_argument("--reference", dest="reference", help="The absolute path to the reference genome.", required=True)
    parser.add_argument("--nanopore", dest="nanopore", help="The path to the nanopore FASTQ.", default=None)
    parser.add_argument("--illumina1", dest="illumina1", help="The path to the first illumina FASTQ.", default=None)
    parser.add_argument("--illumina2", dest="illumina2", help="The path to the second illumina FASTQ.", default=None)
    parser.add_argument("--output", dest="output", help="The path to the output file.", default=None)
    parser.add_argument("--cores", dest="cores", help="The number of cores to use.", default=1)
    args = parser.parse_args()
    if args.nanopore and (args.illumina1 or args.illumina2):
        parser.error("Only one of --nananopore or --illumina1/illumina2 can be specified per run.")
    if args.illumina1 and not args.illumina2:
        parser.error("--illumina2 needs to be specified alongside --illumina1.")
    if args.illumina2 and not args.illumina1:
        parser.error("--illumina1 needs to be specified alongside --illumina2.")
    return args

def map_nanopore_fastq_to_ref(fastq_file,
                            reference_file,
                            output_bam,
                            cores):
    # Build the minimap2 command
    output_sam = output_bam.replace(".bam", ".sam")
    command = f"minimap2 --eqx -t {cores} -a -x map-ont --secondary=no -o {output_sam} {reference_file} {fastq_file} && samtools sort -@ {cores} {output_sam} > {output_bam} && samtools index {output_bam}"
    subprocess.run(command, shell=True, check=True)
    return output_bam

def map_illumina_fastq_to_ref(illumina1_file,
                            illumina2_file,
                            reference_file,
                            output_bam,
                            cores):
    # Build the minimap2 command
    output_sam = output_bam.replace(".bam", ".sam")
    command = f"minimap2 --eqx -t {cores} -a -x sr --secondary=no -o {output_sam} {reference_file} {illumina1_file} {illumina2_file} && samtools sort -@ {cores} {output_sam} > {output_bam} && samtools index {output_bam}"
    subprocess.run(command, shell=True, check=True)
    return output_bam

def get_mean_read_depth_per_contig(bam_file):
    # Run samtools depth command and capture output
    result = subprocess.run(
        ["samtools", "depth", "-aa", bam_file],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=True
    )
    # Read the output into a DataFrame
    data = []
    for line in result.stdout.strip().split("\n"):
        if line != "":
            contig, position, depth = line.split("\t")
            data.append((contig, int(position), int(depth)))
    # Create a DataFrame from the parsed data
    df = pd.DataFrame(data, columns=["contig", "position", "depth"])
    # Calculate mean depth for each contig
    mean_depth_per_contig = df.groupby("contig")["depth"].mean().to_dict()
    return mean_depth_per_contig, data

def count_soft_clipping(bam_file_path):
    # Open the BAM file
    bamfile = pysam.AlignmentFile(bam_file_path, "rb")
    # Dictionary to store soft-clipping counts for each reference position
    clipping_counts = {}
    # Iterate over each read in the BAM file
    for read in bamfile.fetch():
        if read.reference_name not in clipping_counts:
            clipping_counts[read.reference_name] = {}
        if read.cigartuples:
            ref_pos = read.reference_start
            for op, length in read.cigartuples:
                # Soft clipping operation in the CIGAR string is represented by 4
                if op == 4 or op == 5:
                    for i in range(length):
                        if ref_pos + i not in clipping_counts[read.reference_name]:
                            clipping_counts[read.reference_name][ref_pos + i] = 0
                        clipping_counts[read.reference_name][ref_pos + i] += 1
    bamfile.close()
    return clipping_counts

def plot_read_depths(mean_depth_per_contig, data, clipping_counts, coverage_plot):
    # Initialize the plot
    fig, axes = plt.subplots(2, 1, figsize=(20, 10), sharey=False, sharex=True)

    # Plot read depths in the first subplot
    xvals = []
    yvals = []
    current_pos = 0
    prev_contig = None
    v_lines = []
    contig_starts = {}
    sorted_data = sorted(data, key=lambda x: x[0])
    for contig, position, depth in sorted_data:
        if mean_depth_per_contig[contig] > 0:
            xvals.append(current_pos)
            yvals.append(depth / mean_depth_per_contig[contig])
        current_pos += 1
        if contig != prev_contig:
            v_lines.append(current_pos)
            prev_contig = contig
            contig_starts[contig] = current_pos
    axes[0].plot(xvals, yvals, label="Normalized Read Depth", color="blue")
    axes[0].set_ylim([0, 5])
    # Add vertical lines for contig boundaries
    for v in v_lines:
        axes[0].axvline(x=v, linestyle="dashed", color="gray")

    clipping_xvals = []
    clipping_yvals = []
    sorted_contigs = sorted(list(clipping_counts.keys()))
    for contig in sorted_contigs:
        start = contig_starts[contig]
        for position in clipping_counts[contig]:
            clipping_xvals.append(start + position)
            clipping_yvals.append(clipping_counts[contig][position])

    axes[1].scatter(clipping_xvals, clipping_yvals, label="Soft-Clipping Counts", color="red", s=5)

    # Add vertical lines for contig boundaries
    for v in v_lines:
        axes[1].axvline(x=v, linestyle="dashed", color="gray")
    axes[0].grid(True)
    axes[1].grid(True)
    axes[0].set_xlabel("Position")
    axes[1].set_xlabel("Position")
    axes[0].set_ylabel("Normalized Read Depth")
    axes[1].set_ylabel("Soft clipping counts")
    plt.legend()
    plt.xlim(0, 15000)
    # Save the combined plot
    plt.tight_layout()
    plt.savefig(coverage_plot)
    plt.close()

def main():
    # parse the command line arguments
    args = parse_args()
    # make the output dir
    if not os.path.exists(os.path.dirname(args.output)):
        os.mkdir(os.path.dirname(args.output))
    # map the reads to the consensus
    if args.illumina1 and args.illumina2:
        bam_file = map_illumina_fastq_to_ref(
                args.illumina1,
                args.illumina2,
                args.reference,
                args.output.replace(".pdf", ".bam"),
                args.cores)
    if args.nanopore:
        bam_file = map_nanopore_fastq_to_ref(
                args.nanopore,
                args.reference,
                args.output.replace(".pdf", ".bam"),
                args.cores)
    # get the mean depth of each contig
    mean_depth_per_contig, depth_data = get_mean_read_depth_per_contig(bam_file)
    # get soft clipping positions across the genome
    clipping_counts = count_soft_clipping(bam_file)
    # plot the normalised read depths and soft clipping positions
    plot_read_depths(mean_depth_per_contig, depth_data, clipping_counts, args.output)

if __name__ == "__main__":
    main()
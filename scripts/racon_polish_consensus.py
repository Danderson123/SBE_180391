import argparse
import subprocess
import os
import shutil

def parse_args():
    """define args from the command line"""
    parser = argparse.ArgumentParser(description="Map reads to a reference genome and make diagnostic plots to assess assembly quality.")
    parser.add_argument("--reference", dest="reference", help="The absolute path to the reference genome.", required=True)
    parser.add_argument("--nanopore", dest="nanopore", help="The path to the nanopore FASTQ.", default=None)
    parser.add_argument("--illumina1", dest="illumina1", help="The path to the first illumina FASTQ.", default=None)
    parser.add_argument("--illumina2", dest="illumina2", help="The path to the second illumina FASTQ.", default=None)
    parser.add_argument("--output", dest="output", help="The path to the output file.", default=None)
    parser.add_argument("--cores", dest="cores", help="The number of cores to use.", default=1)
    parser.add_argument("--iterations", dest="iterations", type=int, help="The number of racon iterations.", default=1)
    args = parser.parse_args()
    if args.nanopore and (args.illumina1 or args.illumina2):
        parser.error("Only one of --nananopore or --illumina1/illumina2 can be specified per run.")
    if args.illumina1 and not args.illumina2:
        parser.error("--illumina2 needs to be specified alongside --illumina1.")
    if args.illumina2 and not args.illumina1:
        parser.error("--illumina1 needs to be specified alongside --illumina2.")
    return args

def racon_polish(nanopore, sam_file, reference_fasta, output_fasta):
    racon_command = f"racon -u -t 1 --no-trimming -w 500 {nanopore} {sam_file} {reference_fasta} > {output_fasta}"
    subprocess.run(racon_command, shell=True, check=True)

def map_nanopore_reads(output_dir, reference_fasta, nanopore, output_prefix, cores):
    sam_file = os.path.join(output_dir, f"{output_prefix}.mapped.sam")
    bam_file = os.path.join(output_dir, f"{output_prefix}.mapped.bam")

    os.makedirs(output_dir, exist_ok=True)

    map_command = f"minimap2 -a --MD -t {cores} -x map-ont --eqx -o {sam_file} --secondary=no {reference_fasta} {nanopore}"
    subprocess.run(map_command, shell=True, check=True)
    sort_and_index_command = f"samtools sort -@ {cores} {sam_file} > {bam_file} && samtools index {bam_file}"
    subprocess.run(sort_and_index_command, shell=True, check=True)
    return bam_file

def racon_one_iteration_nanopore(nanopore_fastq,
                                output_file,
                                output_dir,
                                reference_fasta,
                                sam_file,
                                sequence_to_polish,
                                polished_sequence,
                                cores):
    # Map reads
    bam_file = map_nanopore_reads(
        output_dir,
        reference_fasta,
        nanopore_fastq,
        sam_file.replace(".mapped.sam", ""),
        cores
    )

    # Polish sequence
    racon_polish(
        nanopore_fastq,
        os.path.join(output_dir, sam_file),
        os.path.join(output_dir, sequence_to_polish),
        os.path.join(output_dir, polished_sequence),
    )

    # Overwrite input for next iteration
    shutil.copy(
        os.path.join(output_dir, polished_sequence),
        os.path.join(output_dir, sequence_to_polish),
    )

    return bam_file

def main():
    args = parse_args()
    # make the output directory if it doesn't exist
    os.makedirs(args.output, exist_ok=True)
    # copy the reference to the output directory
    new_ref = os.path.join(os.path.join(os.path.dirname(args.output), "sequence_to_polish.fasta"))
    shutil.copy(args.reference, new_ref)
    if args.nanopore:
        for i in range(args.iterations):
            racon_one_iteration_nanopore(
                args.nanopore,
                args.output,
                os.path.dirname(args.output),
                new_ref,
                "read.mapped.sam",
                "sequence_to_polish.fasta",
                "polished_sequence.fasta",
                args.cores
            )

if __name__ == "__main__":
    main()
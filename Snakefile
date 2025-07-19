import glob
import os

# list the nanopore reads
nanopore_reads = glob.glob(os.path.join(config["nanopore_directory"], "*.fastq*")) + glob.glob(os.path.join(config["nanopore_directory"], "*.fq*"))
nanopore_samples = [os.path.basename(read).split('.')[0] for read in nanopore_reads]
# list the illumina reads
illumina_reads = glob.glob(os.path.join(config["illumina_directory"], "*.fastq*")) + glob.glob(os.path.join(config["illumina_directory"], "*.fq*"))
illumina_samples = [os.path.basename(read).replace("_1.", ".").replace("_2.", ".").split('.')[0] for read in illumina_reads]

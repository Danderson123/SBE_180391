import glob
import os

# Load the configuration
configfile: 'config.yaml'

# list the nanopore reads
nanopore_reads = glob.glob(os.path.join(config["nanopore_directory"], "*.fastq*")) + glob.glob(os.path.join(config["nanopore_directory"], "*.fq*"))
nanopore_samples = [os.path.basename(read).split('.')[0].replace("_1", "")  for read in nanopore_reads]
# list the illumina reads
illumina_reads = glob.glob(os.path.join(config["illumina_directory"], "*.fastq*")) + glob.glob(os.path.join(config["illumina_directory"], "*.fq*"))
illumina_samples = [os.path.basename(read).replace("_1.", ".").replace("_2.", ".").split('.')[0] for read in illumina_reads]

# define rule all
rule all:
    input:
        expand(os.path.join(config["output_dir"], "{sample}", "samtools", "diagnostic_illumina_plot.pdf"), sample=illumina_samples),
        expand(os.path.join(config["output_dir"], "{sample}", "samtools", "diagnostic_nanopore_plot.pdf"), sample=nanopore_samples),
        expand(os.path.join(config["output_dir"], "{sample}", "racon", "diagnostic_illumina_plot.pdf"), sample=illumina_samples),
        expand(os.path.join(config["output_dir"], "{sample}", "racon", "diagnostic_nanopore_plot.pdf"), sample=nanopore_samples)

# validate the inputs
rule validate_inputs:
    input:
        reference=os.path.join(config["reference_FASTA"]),
        illumina=illumina_reads,
        nanopore=nanopore_reads
    output:
        touch("validation.done")
    run:
        if not os.path.exists(input.reference):
            raise FileNotFoundError(f"Reference file {input.reference} does not exist.")
        # check the nanopore reads
        assert all(r.endswith("_1.fastq.gz") for r in input.nanopore), "All nanopore FASTQ files should end with '_1.fastq.gz'."
        # check the illumina reads
        illumina_samples = set([os.path.basename(read).replace("_1.", ".").replace("_2.", ".").split('.')[0] for read in input.illumina])
        assert all(os.path.exists(os.path.join(config["illumina_directory"], f"{s}_1.fastq.gz")) and os.path.exists(os.path.join(config["illumina_directory"], f"{s}_2.fastq.gz")) for s in illumina_samples), "Illumina reads must be paired and deinterleaved files ending with '_1.fastq.gz' and '_2.fastq.gz'."

# map the Illumina reads to the reference
rule map_illumina_reads:
    input:
        "validation.done",
        reference=os.path.join(config["reference_FASTA"]),
        illumina_1=os.path.join(config["illumina_directory"], "{sample}_1.fastq.gz"),
        illumina_2=os.path.join(config["illumina_directory"], "{sample}_2.fastq.gz")
    output:
        sam_file=os.path.join(config["output_dir"], "{sample}", "mapped_illumina.sam"),
        bam_file=os.path.join(config["output_dir"], "{sample}", "mapped_illumina.bam")
    threads: 4
    shell:
        "minimap2 -ax sr -a --MD --eqx -t {threads} -o {output.sam_file} {input.reference} {input.illumina_1} {input.illumina_2} \
            && samtools sort -@ {threads} {output.sam_file} > {output.bam_file} && samtools index {output.bam_file}"

# map the nanopore reads to the reference
rule map_nanopore_reads:
    input:
        "validation.done",
        reference=os.path.join(config["reference_FASTA"]),
        nanopore=os.path.join(config["nanopore_directory"], "{sample}_1.fastq.gz")
    output:
        sam_file=os.path.join(config["output_dir"], "{sample}", "mapped_nanopore.sam"),
        bam_file=os.path.join(config["output_dir"], "{sample}", "mapped_nanopore.bam")
    threads: 4
    shell:
        "minimap2 -ax map-ont -a --MD --eqx -t {threads} -o {output.sam_file} {input.reference} {input.nanopore} \
            && samtools sort -@ {threads} {output.sam_file} > {output.bam_file} && samtools index {output.bam_file}"

# make illumina consenus using samtools
rule samtools_illumina_consensus:
    input:
        bam_file=os.path.join(config["output_dir"], "{sample}", "mapped_illumina.bam"),
        reference=os.path.join(config["reference_FASTA"])
    output:
        consensus=os.path.join(config["output_dir"], "{sample}", "samtools", "illumina_consensus.fasta")
    threads: 4
    shell:
        "samtools consensus -f FASTA -o {output.consensus} --show-ins no -@ {threads} {input.bam_file}"

# make illumina consenus using racon
rule racon_illumina_consensus:
    input:
        reference=os.path.join(config["reference_FASTA"]),
        illumina_1=os.path.join(config["illumina_directory"], "{sample}_1.fastq.gz"),
        illumina_2=os.path.join(config["illumina_directory"], "{sample}_2.fastq.gz")
    output:
        consensus=os.path.join(config["output_dir"], "{sample}", "racon", "illumina_consensus.fasta")
    threads: 4
    shell:
        "python3 scripts/racon_polish_consensus.py --reference {input.reference} --illumina1 {input.illumina_1} --illumina2 {input.illumina_2} --output {output.consensus} --cores {threads} --iterations 3"

# make nanopore consenus using racon
rule racon_nanopore_consensus:
    input:
        reference=os.path.join(config["reference_FASTA"]),
        nanopore=os.path.join(config["nanopore_directory"], "{sample}_1.fastq.gz")
    output:
        consensus=os.path.join(config["output_dir"], "{sample}", "racon", "nanopore_consensus.fasta")
    threads: 4
    shell:
        "python3 scripts/racon_polish_consensus.py --reference {input.reference} --nanopore {input.nanopore} --output {output.consensus} --cores {threads} --iterations 3"

# make nanopore consenus using samtools
rule samtools_nanopore_consensus:
    input:
        bam_file=os.path.join(config["output_dir"], "{sample}", "mapped_nanopore.bam"),
        reference=os.path.join(config["reference_FASTA"])
    output:
        consensus=os.path.join(config["output_dir"], "{sample}", "samtools", "nanopore_consensus.fasta")
    threads: 4
    shell:
        "samtools consensus -f FASTA -o {output.consensus} --show-ins no -@ {threads} {input.bam_file}"

# make diagnostic plots for the illumina reads
rule diagnostic_samtools_illumina:
    input:
        consensus=os.path.join(config["output_dir"], "{sample}", "samtools", "illumina_consensus.fasta"),
        illumina_1=os.path.join(config["illumina_directory"], "{sample}_1.fastq.gz"),
        illumina_2=os.path.join(config["illumina_directory"], "{sample}_2.fastq.gz")
    output:
        plot=os.path.join(config["output_dir"], "{sample}", "samtools", "diagnostic_illumina_plot.pdf")
    threads: 4
    shell:
        "python3 scripts/make_diagnostic_plots.py --reference {input.consensus} --illumina1 {input.illumina_1} --illumina2 {input.illumina_2} --output {output.plot} --cores {threads}"

# make diagnostic plots for the nanopore readsillumina_consensus
rule diagnostic_samtools_nanopore:
    input:
        consensus=os.path.join(config["output_dir"], "{sample}", "samtools", "nanopore_consensus.fasta"),
        nanopore=os.path.join(config["nanopore_directory"], "{sample}_1.fastq.gz")
    output:
        plot=os.path.join(config["output_dir"], "{sample}", "samtools", "diagnostic_nanopore_plot.pdf")
    threads: 4
    shell:
        "python3 scripts/make_diagnostic_plots.py --reference {input.consensus} --nanopore {input.nanopore} --output {output.plot} --cores {threads}"

# make diagnostic plots for the illumina reads
rule diagnostic_racon_illumina:
    input:
        consensus=os.path.join(config["output_dir"], "{sample}", "racon", "illumina_consensus.fasta"),
        illumina_1=os.path.join(config["illumina_directory"], "{sample}_1.fastq.gz"),
        illumina_2=os.path.join(config["illumina_directory"], "{sample}_2.fastq.gz")
    output:
        plot=os.path.join(config["output_dir"], "{sample}", "racon", "diagnostic_illumina_plot.pdf")
    threads: 4
    shell:
        "python3 scripts/make_diagnostic_plots.py --reference {input.consensus} --illumina1 {input.illumina_1} --illumina2 {input.illumina_2} --output {output.plot} --cores {threads}"

# make diagnostic plots for the nanopore readsillumina_consensus
rule diagnostic_racon_nanopore:
    input:
        consensus=os.path.join(config["output_dir"], "{sample}", "racon", "nanopore_consensus.fasta"),
        nanopore=os.path.join(config["nanopore_directory"], "{sample}_1.fastq.gz")
    output:
        plot=os.path.join(config["output_dir"], "{sample}", "racon", "diagnostic_nanopore_plot.pdf")
    threads: 4
    shell:
        "python3 scripts/make_diagnostic_plots.py --reference {input.consensus} --nanopore {input.nanopore} --output {output.plot} --cores {threads}"

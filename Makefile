# Makefile for SBE_180391 project
PROJECT = SBE_180391
OS := $(shell uname -s)
CORES = 16

.PHONY: run snakemake-run singularity-run test

# Standard snakemake execution
run:
	snakemake --cores $(CORES) --nolock --rerun-incomplete

# Snakemake with Singularity container
singularity:
	singularity run $(PROJECT).img snakemake --cores $(CORES) --nolock --rerun-incomplete

# Dry run for testing
test:
	snakemake --dry-run --cores $(CORES)

install:
	conda env create -f environment.yaml
	conda activate $(PROJECT)

# build Singularity image
singularity-build:
	sudo singularity build $(PROJECT).img Singularity.def
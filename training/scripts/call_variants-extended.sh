#!/usr/bin/env bash

# set bash strict mode
set -euo pipefail

# allow debug mode by running `TRACE=1 ./map.sh`
if [[ "${TRACE-0}" == "1" ]]; then set -x; fi

######################################################################
# Script to call variants and create vcf files from sorted bam files #
######################################################################

# get file path of script
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
# Otherwise, you would need to make sure to call the script from within the
# directory where it is stored.
# Alternatively, use absolute paths, but this makes the script less portable.

# set number of threads for downstream tools
n_threads=2

# define in and outputs and create directories when necessary
# by defining everything here, there is no need to `cd` to directories first
reference_genome="${SCRIPT_DIR}/../data/reference/PlasmoDB-67_Pfalciparum3D7_Genome.fasta"
alignment_dir="${SCRIPT_DIR}/../results/dups"
output_dir="${SCRIPT_DIR}/../results/gatk"
mkdir -p "${output_dir}"

# create fai and dict files
samtools faidx "${reference_genome}"
if ! [ -f "${reference_genome%.fasta}.dict" ]; then
	echo "${reference_genome%.fasta}.dict"
	gatk CreateSequenceDictionary -R "${reference_genome}"
fi

# loop through bam files and call variants using HaplotypeCaller
for bam in "${alignment_dir}"/*.bam; do
	sample_name=$(basename ${bam} .bam)

	gatk HaplotypeCaller \
		-R "${reference_genome}" \
		-I "${bam}" \
		-O "${output_dir}/${sample_name}.g.vcf.gz" \
		--native-pair-hmm-threads ${n_threads} \
		-ERC GVCF
done

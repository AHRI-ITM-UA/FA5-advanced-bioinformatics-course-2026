#!/usr/bin/env bash

# set bash strict mode
set -euo pipefail

# allow debug mode by running `TRACE=1 ./map.sh`
if [[ "${TRACE-0}" == "1" ]]; then set -x; fi

##########################################################
# Script to map fastq files to reference genome with bwa #
##########################################################

# get file path of script
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
# Otherwise, you would need to make sure to call the script from within the
# directory where it is stored.
# Alternatively, use absolute paths, but this makes the script less portable.

# set number of threads for downstream tools
n_threads=2

# define in and outputs and create directories when necessary
# by defining everything here, there is no need to `cd` to directories first
ref="${SCRIPT_DIR}/../data/reference/PlasmoDB-67_Pfalciparum3D7_Genome.fasta"
fastq_dir="${SCRIPT_DIR}/../results/trimmomatic"
output_dir="${SCRIPT_DIR}/../results/bwa/"
mkdir -p "${output_dir}"

# check if fastq directory exist
if [ ! -d "${fastq_dir}" ]; then
    echo "FASTQ input directory (${fastq_dir}) does not exist."
fi

# check if reference fasta file exists
if ! [ -f "${ref}" ]; then
    echo "Reference fasta file not found (${ref})."
fi

# log run options
printf "
BWA MEM script | $(basename "$0")
==============================================

Output directory:           ${output_dir}
FASTQ reads directory:      ${fastq_dir}
Reference pfalciparum:      ${ref}
threads:                    ${n_threads}
"

###################
# start of script #
###################

# create a bwa index file
bwa index "${ref}"

# loop through read pairs and map them using bwa
for read_1 in "${fastq_dir}"/*_R1_001.fastq.gz; do

    # get sample name for each read pair
    sample_name=$(basename ${read_1} _R1_001.fastq.gz)

    printf "\nProcessing sample %s ...\n\n" "${sample_name}"

    # map sample reads to reference genome
    bwa mem \
        -t "${n_threads}" \
        -Y -K 100000000 \
        -R "@RG\tID:L001\tPL:ILLUMINA\tSM:${sample_name}" \
        "${ref}" \
        "${fastq_dir}/${sample_name}_R1_001.fastq.gz" \
        "${fastq_dir}/${sample_name}_R2_001.fastq.gz" |
        # sort and compress to bam
        samtools sort \
            -@ "${n_threads}" \
            -o "${output_dir}/${sample_name}.sort.bam"

    # generate samtools report
    samtools flagstat \
        -@ ${n_threads} \
        "${output_dir}/${sample_name}.sort.bam" \
        >"${output_dir}/${sample_name}.flagstat"

    # create bam index .bai
    samtools index \
        -@ ${n_threads} \
        "${output_dir}/${sample_name}.sort.bam"
done

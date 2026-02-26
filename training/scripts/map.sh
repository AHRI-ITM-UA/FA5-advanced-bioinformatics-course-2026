#!/usr/bin/env bash

##########################################################
# Script to map fastq files to reference genome with bwa #
##########################################################

# make sure to run this script from within the directory where it is stored!

# move to the directory containing the fastq files
cd ../results/trimmomatic

# create a directory to store the results and store the path as a variable
output_dir="../bwa-simple"
mkdir -p ${output_dir}

# create a bwa index file
bwa index ../../data/reference/PlasmoDB-67_Pfalciparum3D7_Genome.fasta

# loop through read pairs and map them using bwa
for read_1 in *_R1_001.fastq.gz; do
    sample_name=$(basename ${read_1} _R1_001.fastq.gz)

    # map sample reads to reference genome
    bwa mem \
        -t 2 \
        -Y -K 100000000 \
        -R "@RG\tID:L001\tPL:ILLUMINA\tSM:${sample_name}" \
        ../../data/reference/PlasmoDB-67_Pfalciparum3D7_Genome.fasta \
        ${read_1} \
        ${sample_name}_R2_001.fastq.gz |
        samtools sort -@ 4 -o ${output_dir}/${sample_name}.sort.bam

    # generate samtools report
    samtools flagstat -@ 4 ${output_dir}/${sample_name}.sort.bam \
        >${output_dir}/${sample_name}.flagstat

    # create bam index .bai
    samtools index \
        -@ 2 \
        "${output_dir}/${sample_name}.sort.bam"

done

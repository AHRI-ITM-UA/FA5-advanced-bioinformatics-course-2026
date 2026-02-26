#!/usr/bin/env bash

######################################################################
# Script to call variants and create vcf files from sorted bam files #
######################################################################

# make sure to run this script from within the directory where it is stored!

# store path to reference genome
reference_genome="../data/reference/PlasmoDB-67_Pfalciparum3D7_Genome.fasta"

# create a directory to store the results and store the path as a variable
output_dir="../results/gatk"
mkdir -p ${output_dir}

# create fai and dict files
samtools faidx ${reference_genome}
gatk CreateSequenceDictionary -R ${reference_genome}

# loop through bam files and call variants using gatk
for bam in ../results/dups/*.bam; do
	sample_name=$(basename ${bam} .bam)
	gatk HaplotypeCaller \
		-R ${reference_genome} \
		-I ${bam} \
		-O ${output_dir}/${sample_name}.g.vcf.gz \
		--native-pair-hmm-threads 2 \
		-ERC GVCF
done

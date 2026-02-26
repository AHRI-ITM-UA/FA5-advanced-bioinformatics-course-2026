## activate the conda environment
conda activate FA5_ethiopia_2026

## vcf2phylip not working on macbook
export PATH=/Users/pmonsieurs/programming/software/vcf2phylip:/Users/pmonsieurs/programming/software/plink2:/Users/pmonsieurs/programming/software/admixture/:${PATH}

############ 1. Background ##############

############ 2. Explore the input data ##############

## specifying some input output directories. Those can be changed to match your own directory structure. You are not obliged to use those parameters, so you can also "hard-code" everyhting Also specify the number of threads to use for the analyses, which can speed up the analyses if you have a multi-core processor.

data_dir=/Users/pmonsieurs/programming/ext_FA5_ethiopia_2026/data/
ref_genome_dir=${data_dir}/ref_genome/
ref_genome=${ref_genome_dir}/Ldonovani_BPK282.fasta
results_dir=/Users/pmonsieurs/programming/ext_FA5_ethiopia_2026/results/
threads=4


## download the data and check the amount of data in the fastq file
cd ${data_dir}/fastq/
for file in *.fastq.gz; do \
    echo $file\'
    wc -l $file
done


############ 3. Quality control  ##############

## start with running fastQC
mkdir -p ${results_dir}/fastqc/
fastqc ${data_dir}/fastq/*.fastq.gz -o ${results_dir}/fastqc/

## combine individual fastQC reports into a single report using multiQC
multiqc ${results_dir}/fastqc/ -o ${results_dir}/fastqc/

## optional: if your data contain adapters in the multiQC report, you can run fastp to trim the adapters from the raw data. But most modern alignment algorithms run well even with adapters in the data, as they do soft clipping, which means that they will ignore the adapter sequences during the alignment process (and also low-quality bases at the ends of the reads). 
output_dir=${data_dir}/fastq/fastq_trimmed/
mkdir -p ${output_dir}
threads=4

for fastq_file_R1 in ${data_dir}/fastq/*_R1.fastq.gz; do
    
    # Get base name (remove path and _R1.fastq.gz)
    base=$(basename ${fastq_file_R1} _R1.fastq.gz)
    
    fastq_file_R2=${data_dir}/fastq/${base}_R2.fastq.gz

    echo "Processing ${base}..."

    fastp \
        -i ${fastq_file_R1} \
        -I ${fastq_file_R2} \
        -o ${output_dir}/${base}_R1.trimmed.fastq.gz \
        -O ${output_dir}/${base}_R2.trimmed.fastq.gz \
        --detect_adapter_for_pe \
        --thread ${threads} \
        --html ${output_dir}/${base}.fastp.html \
        --json ${output_dir}/${base}.fastp.json

done



########### 4. MAPPING WITH BWA-MEM AND SAMTOOLS #############

## run BWA to align the reads to the reference genome. Also set some
## parameters for BWA, such as the seed length and the number of threads to use.
seed=30
mapq_cutoff=30
mkdir -p ${results_dir}/bwa/

## index the reference genome using BWA, which is needed for the alignment step. This will create several index files that are used by BWA to efficiently align the reads to the reference genome. The indexing step only needs to be done once for each reference genome, and the resulting index files can be reused for multiple alignments.
bwa index ${ref_genome}


# for fastq_file_R1 in ${data_dir}/fastq/*_R1.fastq.gz; do
for fastq_file_R1 in ${data_dir}/fastq/fastq_trimmed/*_R1.trimmed.fastq.gz; do
    
    ## if you do without trimming
    # sample_name=$(basename ${fastq_file_R1} _R1.fastq.gz)
    # fastq_file_R2=${fastq_file_R1/_R1.fastq.gz/_R2.fastq.gz}

    ## if you do with trimming using fastp
    sample_name=$(basename ${fastq_file_R1} _R1.trimmed.fastq.gz)
    fastq_file_R2=${fastq_file_R1/_R1.trimmed.fastq.gz/_R2.trimmed.fastq.gz}

    #### 1) Align reads to the reference genome using BWA-MEM and \
    #### sort the output BAM file using samtools sort. Also add read group information to the BAM file using the -R option of BWA-MEM.
    bwa mem \
        -R "@RG\tID:${sample_name}\tSM:${sample_name}\tPL:ILLUMINA" \
        -k $seed \
        -t $threads \
        ${ref_genome} \
        ${fastq_file_R1} \
        ${fastq_file_R2} | \
    samtools sort -@$threads -o ${results_dir}/bwa/${sample_name}.bam
    samtools flagstat ${results_dir}/bwa/${sample_name}.bam > ${results_dir}/bwa/${sample_name}.flagstat

    ### 2) Filter on mapping quality using samtools view, and keep only reads with a mapping quality above a certain threshold (e.g., 30). This is important to ensure that the downstream analyses are based on high-quality alignments, and to avoid biases in the results due to low-quality reads that may be misaligned or have multiple mapping locations in the genome.
    samtools view -@ $threads -b -q $mapq_cutoff ${results_dir}/bwa/${sample_name}.bam -o ${results_dir}/bwa/${sample_name}.mapq${mapq_cutoff}.bam
    samtools index -@ $threads ${results_dir}/bwa/${sample_name}.mapq${mapq_cutoff}.bam
    samtools flagstat ${results_dir}/bwa/${sample_name}.mapq${mapq_cutoff}.bam > ${results_dir}/bwa/${sample_name}.mapq${mapq_cutoff}.flagstat


    ### 3) Remove duplicates using Picard MarkDuplicates, and keep only properly paired reads using samtools view with the -f 0x2 option. Removing duplicates is important to avoid biases in the downstream analyses due to PCR amplification artifacts, which can lead to overrepresentation of certain reads and affect the accuracy of variant calling and other analyses. 
    picard MarkDuplicates \
        REMOVE_DUPLICATES=true \
        I=${results_dir}/bwa/${sample_name}.mapq${mapq_cutoff}.bam \
        O=${results_dir}/bwa/${sample_name}.mapq${mapq_cutoff}.rmdup.bam \
        M=${results_dir}/bwa/${sample_name}.mapq${mapq_cutoff}.rmdup.metrics.txt

    samtools index -@ $threads ${results_dir}/bwa/${sample_name}.mapq${mapq_cutoff}.rmdup.bam
    samtools flagstat ${results_dir}/bwa/${sample_name}.mapq${mapq_cutoff}.rmdup.bam > ${results_dir}/bwa/${sample_name}.mapq${mapq_cutoff}.rmdup.flagstat


    ### 4) Keep only properly paired reads using samtools view with the -f 0x2 option. Properly paired reads are those where both reads in a pair are mapped to the reference genome in the expected orientation and distance from each other. Keeping only properly paired reads can help to improve the accuracy of variant calling and other analyses, as these reads are more likely to represent true biological signals rather than artifacts or noise in the data.
    samtools view -@ $threads -b -f 0x2 ${results_dir}/bwa/${sample_name}.mapq${mapq_cutoff}.rmdup.bam -o ${results_dir}/bwa/${sample_name}.mapq${mapq_cutoff}.rmdup.proper.bam
    samtools index -@ $threads ${results_dir}/bwa/${sample_name}.mapq${mapq_cutoff}.rmdup.proper.bam
    samtools flagstat ${results_dir}/bwa/${sample_name}.mapq${mapq_cutoff}.rmdup.proper.bam > ${results_dir}/bwa/${sample_name}.mapq${mapq_cutoff}.rmdup.proper.flagstat

done


############# 5. SNP prediction using GATK ##############
gatk CreateSequenceDictionary -R ${ref_genome}
samtools faidx ${ref_genome}

mkdir -p ${results_dir}/gatk/

## 1) Run GATK HaplotypeCaller to call variants for each sample separately, and output the results in GVCF format. This allows us to combine the results from all samples later on.
for bam_file in ${results_dir}/bwa/*.mapq${mapq_cutoff}.rmdup.proper.bam; do
    sample_name=$(basename ${bam_file} .mapq${mapq_cutoff}.rmdup.proper.bam)
    gatk HaplotypeCaller \
        -R ${ref_genome} \
        -I ${bam_file} \
        -O ${results_dir}/gatk/${sample_name}.g.vcf.gz \
        -ERC GVCF
done





## 2) Combine the GVCF files from all samples using GATK CombineGVCFs, and then run GenotypeGVCFs to perform joint genotyping and output a single VCF file with the variants called across all samples.


####  2.1) genomicsDBImport ####
## Option 1: use genomicsDBImport to import the GVCF files into a GenomicsDB, and then run GenotypeGVCFs on the GenomicsDB. This can be more efficient than combining the GVCF files into a single file, especially for large datasets, as it allows for parallel processing and reduces the memory requirements of the analysis. This is also the recommended approach by GATK for large datasets, as it can improve the speed and scalability of the analysis.

## create a sample map file that lists the sample names and the corresponding GVCF files, which is needed for the GenomicsDBImport step. The sample map file should have two columns: the first column contains the sample names, and the second column contains the paths to the corresponding GVCF files.

## option 1: create the sample map using a fro loop
for f in ${results_dir}/gatk/*.g.vcf.gz; do
    echo -e "$(bcftools query -l "$f")\t$f"
done > ${results_dir}/gatk/sample_map.txt

## option 2: create file with sample names and file with gvcf files and then paste them together. Only works if only one sample per gvcf file, and if the sample names are in the same order as the gvcf files. --> more tricky
ls ${results_dir}/gatk/*.g.vcf.gz > ${results_dir}/gatk/files.txt
for file in ${results_dir}/gatk/*.g.vcf.gz; do
    bcftools query -l "$file"
done > ${results_dir}/gatk/samples.txt
paste ${results_dir}/gatk/samples.txt ${results_dir}/gatk/files.txt > ${results_dir}/gatk/sample_map_tricky.txt

## double check that no duplicates exist in the sample map file, and that all samples are included. If duplicates exist, genomicsDBImport will fail, as it requires unique sample names to import the GVCF files into the GenomicsDB. 
cut -f1 ${results_dir}/gatk/sample_map.txt | sort | uniq -d


## run the genomicsDBImport step to import the GVCF files into a GenomicsDB. This will create a GenomicsDB workspace that contains the variant data from all samples, which can be used for joint genotyping in the next step.
grep "^>" ${ref_genome} | sed 's/>//' > ${results_dir}/gatk/genome.intervals
gatk GenomicsDBImport \
    --reference ${ref_genome} \
    --genomicsdb-workspace-path ${results_dir}/gatk/combined \
    --sample-name-map ${results_dir}/gatk/sample_map.txt \
    --reader-threads $threads \
    --intervals ${results_dir}/gatk/genome.intervals

## run the GenotypeGVCFs step on the genomicsDB database
gatk GenotypeGVCFs \
    -R ${ref_genome} \
    -V gendb://${results_dir}/gatk/combined \
    -O ${results_dir}/gatk/combined.vcf.gz


## option 2: use the --variant flag of CombineGVCFs to specify the GVCF files to combine. This is easier to set up than the GenomicsDBImport approach, but it may not be as efficient for large datasets, as it requires combining all the GVCF files into a single file before running GenotypeGVCFs. This approach can be suitable for smaller datasets.
gatk CombineGVCFs \
    -R ${ref_genome} \
    $(for f in ${results_dir}/gatk/*.g.vcf.gz; do echo --variant $f; done) \
    -O ${results_dir}/gatk/combined.old_way.g.vcf.gz


gatk GenotypeGVCFs \
    -R ${ref_genome} \
    -V ${results_dir}/gatk/combined.old_way.g.vcf.gz \
    -O ${results_dir}/gatk/combined.old_way.vcf.gz


## 3) Filter the variants based on quality and other criteria using bcftools filter. For example, we can filter out variants with low quality scores or those that are not SNPs. Also, only keep the biallelic SNPs (i.e., those with only two alleles) and with a quality score above 30. This is needed when you want to go for population genetic analyses, as these analyses typically require high-quality SNPs that are likely to be true positives, and which are only biallelic (i.e., have only two alleles) to avoid complications in the analyses.
bcftools view \
    -m2 -M2 \
    -i 'TYPE="snp" && QUAL>30' \
    ${results_dir}/gatk/combined.vcf.gz \
    -Oz -o ${results_dir}/gatk/combined.filtered.vcf.gz

## count the number of SNPs before and after filtering, to see how many SNPs were retained for downstream analyses. This can be done using bcftools stats or by simply counting the number of lines in the VCF file (excluding the header lines).
bcftools stats ${results_dir}/gatk/combined.vcf.gz | grep "^SN"
bcftools stats ${results_dir}/gatk/combined.filtered.vcf.gz | grep "^SN"

gunzip -c  ${results_dir}/gatk/combined.vcf.gz | grep -v "^#" | wc
gunzip -c  ${results_dir}/gatk/combined.filtered.vcf.gz | grep -v "^#" | wc






############# 4. Phylogenetic analysis ##############

## 1) convert the vcf file to phylip or fasta
mkdir -p ${results_dir}/plink/

cd ${results_dir}/plink/
vcf2phylip.py \
    -i ${results_dir}/gatk/combined.filtered.vcf.gz \
    --fasta


## 2) run IQ-TREE to build a phylogenetic tree based on the SNP data. The model of evolution to use, can be selected using the built-in model selection in IQ-TREE, which will test different models and select the best one based on the data. Also perform bootstrapping to assess the support for the branches in the tree.
iqtree \
    -s ${results_dir}/plink/combined.filtered.min4.fasta \
    -m MFP \
    -bb 1000 \
    -nt AUTO \
    -o Lmaj_Friedlin

## optional: remove the outgroup and rerun the threee construction to see how the tree topology changes. This can be done by simply removing the outgroup sequence from the fasta file and then rerunning the IQ-TREE analysis on the modified fasta file. Another option could be to delete the sample using bcftools view -s ^Lmaj_Friedlin
cp combined.filtered.min4.fasta combined.filtered.min4.no_outgroup.fasta

iqtree \
    -s ${results_dir}/plink/combined.filtered.min4.no_oufasta \
    -m MFP \
    -bb 1000 \
    -nt AUTO 

## 3) alternative: run RAxML-NG to build a phylogenetic tree based on the SNP data. The model of evolution to use, you have to select yourself. Also perform bootstrapping to assess the support for the branches in the tree.
raxml-ng --all \
  --msa ${results_dir}/plink/combined.filtered.min4.fasta \
  --model GTR+G \
  --bs-trees 10 \
  --threads auto

## 4) alternative 2: run FastTree to build a phylogenetic tree based on the SNP data. The model of evolution to use, you have to select yourself. FastTree is a very fast method for building phylogenetic trees, but it may not be as accurate as IQ-TREE or RAxML-NG, especially for large datasets or complex models of evolution.
FastTree -gtr -nt ${results_dir}/plink/combined.filtered.min4.fasta > ${results_dir}/plink/combined.filtered.FastTree.tree

## 5) alternative 4: read in fasta into SplitsTree and build a phylogenetic network based on the SNP data. This can be useful to visualize the relationships between the samples, especially if there is evidence of recombination or gene flow between them, which can complicate the interpretation of a traditional phylogenetic tree. The model of evolution to use, you have to select yourself.



########### 5. PCA analysis ##############
mkdir -p ${results_dir}/plink/
cp ${data_dir}/vcf_full/full.filtered.vcf.gz ${results_dir}/gatk/full.filtered.vcf.gz

## run PLINK to convert the VCF file to PLINK format, and then perform some quality control on the SNPs and individuals. For example, we can filter out SNPs that are missing in more than 5% of the samples, or those with a minor allele frequency (MAF) below 5%, and we can also filter out individuals that are missing more than 10% of the SNPs. This is important to ensure that the PCA analysis is based on high-quality data and to avoid biases in the results.

plink2 --vcf ${results_dir}/gatk/full.filtered.vcf.gz --make-bed --out ${results_dir}/plink/full.filtered
# plink2 --vcf ${results_dir}/gatk/full.filtered.withIDs.vcf.gz --make-bed --out ${results_dir}/plink/full.filtered

# geno filters on SNPs: --geno 0.05 (remove SNPs with >5% missing data), --maf 0.05 (remove SNPs with MAF < 5%)
# mind filter on individuals: --mind 0.1 (remove individuals with >10% missing data)
plink2 --bfile ${results_dir}/plink/full.filtered \
      --geno 0.05 \
      --maf 0.05 \
      --mind 0.1 \
      --make-bed \
      --out ${results_dir}/plink/full.filtered.QC

## include LD pruning to remove SNPs that are in high linkage disequilibrium (LD) with each other, which can bias the PCA results. This is typically done by calculating the pairwise LD between SNPs and removing one of the SNPs in each pair that has an LD above a certain threshold (e.g., r^2 > 0.2). After LD pruning, we can perform PCA on the pruned dataset to visualize the genetic structure of the samples and identify any clusters or patterns in the data.
# plink2 \
#   --bfile ${results_dir}/plink/full.filtered.QC \
#   --set-all-var-ids @:# \
#   --make-bed \
#   --out ${results_dir}/plink/full.filtered.QC.uniqueIDs

plink2 \
  --bfile ${results_dir}/plink/full.filtered.QC \
  --indep-pairwise 50 5 0.5 \
  --out ${results_dir}/plink/full.filtered.QC.pruned \
  --set-all-var-ids \@:# \
  --bad-ld

plink2 \
  --bfile ${results_dir}/plink/full.filtered \
  --extract ${results_dir}/plink/full.filtered.QC.pruned.prune.in \
  --make-bed \
  --set-all-var-ids \@:# \
  --out ${results_dir}/plink/full.filtered.QC.pruned

plink2 \
  --bfile ${results_dir}/plink/full.filtered.QC.pruned \
  --freq \
  --out ${results_dir}/plink/full.filtered.QC.pruned

## run PCA. PCA allows to visualize the genetic structure of the samples and identify any clusters or patterns in the data. The PCA results can be visualized using a scatter plot, where each point represents a sample and the axes represent the principal components. Samples that are genetically similar will cluster together in the plot, while those that are genetically different will be more spread out. This can help to identify any population structure or substructure in the data, as well as any outliers or unusual samples that may require further investigation. (limit to 10 principal components)
plink2 \
  --bfile ${results_dir}/plink/full.filtered.QC.pruned \
  --read-freq ${results_dir}/plink/full.filtered.QC.pruned.afreq \
  --pca 5 \
  --out ${results_dir}/plink/full.filtered.QC.pruned.PCA

## do visulisation in R --> see phylo_viz.R


########### 6. Admixture analysis ##############

## run admixture to estimate the ancestry proportions of each sample, which can provide insights into the population structure and history of the samples. Admixture is a model-based clustering method that estimates the proportion of ancestry from different ancestral populations for each sample. The number of ancestral populations (K) can be specified by the user, and the algorithm will assign each sample to one or more of these populations based on their genetic data. The results can be visualized using a bar plot, where each bar represents a sample and the colors represent the different ancestral populations. This can help to identify any admixture events or gene flow between populations, as well as any patterns of ancestry that may be relevant for understanding the evolutionary history of the samples.

## convert bim file
awk '{$1="0";print $0}' ${results_dir}/plink/full.filtered.QC.bim > ${results_dir}/plink/full.filtered.QC.bim.tmp
mv ${results_dir}/plink/full.filtered.QC.bim ${results_dir}/plink/full.filtered.QC.bim.source
mv ${results_dir}/plink/full.filtered.QC.bim.tmp ${results_dir}/plink/full.filtered.QC.bim

## run for multiple values of K
for K in {2..10}; do
    echo "Running admixture for K=${K}..."
    admixture --cv ${results_dir}/plink/full.filtered.QC.bed $K > admixture_log${K}.out
done

## select the number of populations (K) with the lowest cross-validation error
grep "CV error" ${results_dir}/admixture_log*.out

## do visulisation in R --> see phylo_viz.R


########### 7. IBD analysis ############

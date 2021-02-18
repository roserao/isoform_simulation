Step 1: Randomly Select Genes

Source Code:
get_gene_annotation.R
Sample 10 genes for each number of isoforms

	Input: 
	1_gene_list_filtered.csv: all the isoforms of the genes in the list are
	expressed in our data. That is, no NA value will exist in the isoform 
	expression covariance matrix.

	Output:
	gene_count.csv
	under the directory gene_sample
		directory {1-10}: each gene has a csv file recording info of its
		isoforms
		gene{1-10}.csv: "gene_id", "start", "end", "chr"

Upload gene{1-10}.csv to Hoffman cluster

Step 2: Genotype Simulation

Source Code:
geno_sim.sh:
Extract vcf file of genotype of African population from 1000-Genome Database
Calculate LD matrix using Plink
Simulate genotype for 1000 individuals (2_geno_sim.py)

	Input:
	gene{1-10}.csv: "gene_id", "start", "end", "chr"

	Output:
	./geno/{1-10}/{gene_name} contains
	1) {gene_name}_raw.vcf
	   {gene_name}_AFR.recode.vcf #extract African population
	   {gene_name}_AFR.recode.ped + map #plink format
	2) {gene_name}_AFR.clean.bed + fam + bim #maf filter > 0.05
	3) {gene_name}_AFR.clean.ld # LD matrix V
	4) {gene_name}_AFR_sample.txt # simulated genotype for 1000 individuals

geno_sim.py:
Simulate genotype
Argument:
--gene: name of gene
--size: number of individual genotype generated 
--dir: input and output directory

Step 3: Extract Isoform Expression Data

Source Code:
gene_exp_cov.sh:
Submit job for 3_gene_exp_cov.R

gene_exp_cov.R:
For each gene with a folder, output tpm file of isoform expression

	Output:
	{gene_name}_tpm.csv

Step 4: Isoform Expression Simulation

Source Code:
isoexp_sim.sh
For each gene with a folder, simulate isoform expression with assigned gene expression heritability 0.1 0.2 0.3 0.4 0.5

model.py:
Simulate Isoform Expression
Pipeline described in the report



























 	

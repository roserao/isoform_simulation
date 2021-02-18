library("tximport")
library("rhdf5")
library("dplyr")

common_path <- "/u/project/zarlab/yiwen99/kallisto_quant/results"
samples <- list.files(common_path)
files <- file.path(common_path, samples, "abundance.h5")
names(files) <- paste0(samples)

to_gene <- read.delim(file.path("/u/project/zarlab/hjp/geuvadis_data/annotation",
                                "transcripts_to_genes.txt"), 
                      header = FALSE, stringsAsFactors = FALSE)
colnames(to_gene) <- c("transcript_id", "ens_gene", "common_gene")
tx2gene <- to_gene[, c(1, 2)]
to_gene <- mutate(group_by(to_gene, ens_gene), n_transcripts = length(transcript_id))

# report simply counts, write to a rds file
kal <- tximport(files, type = "kallisto", tx2gene = tx2gene, txOut = TRUE)
tpm <- kal[]

ddir <- "/u/scratch/r/roserao/data/gene_sample/gene"
sdir <- "/u/scratch/r/roserao/data/geno"
for (ntrans in 2:10) {
  print("##################")
  print(ntrans)

  gene_list <- read.csv(paste(ddir, ntrans, ".csv", sep = ""))
  cur_subset = filter(to_gene, n_transcripts == ntrans)

  for (gene in gene_list$gene_id) {
    
    gdir <- file.path(sdir, ntrans, gene)
    if (!dir.exists(gdir)) next

    print(gene)

    cur_transcripts = filter(cur_subset, grepl(gene, ens_gene))$transcript_id
    cur_tpm = kal$abundance[cur_transcripts, ]
    cor_m <- cor(t(cur_tpm))
    print(cor_m)

    if(sum(is.na(cor_m)) > 0) {
      print("NA value exists")
    } else {
      e_check = all(eigen(cor_m)$values > 0)
      print("eigenvalue checking:")
      print(e_check)
    }
    
    #write.csv(cor_m, file.path(gdir, paste(gene, "_isoform_cov.csv", sep = "")))
    write.csv(cur_tpm, file.path(gdir, paste(gene, "_tpm.csv", sep = "")))
  }
  
}

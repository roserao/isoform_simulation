library("biomaRt")

##### pull out gene and isoform data from biomaRt
ensembl <- useMart("ensembl")
# listEnsembl(GRCh=37)
grch37 <- useEnsembl(biomart = "ensembl", GRCh = 37)
# listDatasets(grch37)
ds <- useDataset("hsapiens_gene_ensembl", mart = grch37)
gene_info <- getBM(
  c("ensembl_gene_id", "chromosome_name", "start_position", "end_position",
    "strand", "transcription_start_site", "external_gene_name", "transcript_count"),
  mart = ds
)
# summary(gene_info)

##### filtered gene list with
gene_list_filtered <- read.csv("1_gene_list_filtered.csv")
gene_list_filtered$gene_name <- substr(gene_list_filtered$gene_name, 1, 15)
# table(gene_list_filtered$num_trans)

##### sample 10 genes for each isoform
valid_chro <- as.character(1:22) # we don't need genes on X and Y
region <- 100000
# sample 10 genes for each number of isoforms
for (ntrans in 2:10) {
  gene_sub <- dplyr::distinct(dplyr::filter(gene_info, 
                                            transcript_count == ntrans,
                                            chromosome_name %in% valid_chro,
                                            ensembl_gene_id %in% gene_list_filtered$gene_name))

  dir.create(paste("./gene_sample/", ntrans, sep = ""))

  set.seed(105104394)
  gene_name <- sample(unique(gene_sub$ensembl_gene_id), size = 19)

  i <- 0
  list <- NULL
  for (gene in gene_name) {
    
    if (i == 10) break # only need 10 samples
    
    info <- dplyr::distinct(dplyr::filter(gene_sub, ensembl_gene_id == gene))
    if (nrow(info) == ntrans) {
      
      write.csv(info, paste("./gene_sample/", ntrans, "/", gene, ".csv", sep = ""),
                row.names = FALSE)
      df <- data.frame(
        gene_id = gene,
        start = min(info$transcription_start_site) - region,
        end = max(info$transcription_start_site) + region,
        chr = info$chromosome_name[1]
      )
      list <- rbind(list, df)
      i <- i + 1
    } 
  }
  write.csv(list, paste("./gene_sample/gene", ntrans, ".csv", sep = ""),
    row.names = FALSE
  )
}

##### create summary csv for future notation
list_all <- NULL
for (ntrans in 2:10) {
  org <- read.csv(paste("./gene_sample/gene", ntrans, ".csv", sep = ""))
  org$ntrans <- ntrans
  list_all <- rbind(list_all, org)
}
write.csv(list_all, "gene_count.csv", row.names = FALSE)

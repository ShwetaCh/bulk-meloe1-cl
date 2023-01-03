library(tximportData)
library(data.table)
library(tximport)
library(rhdf5)
library(ensembldb)
library(factoextra)
library(dplyr)
library(DESeq2)
library(limma)
library(ggplot2)

run_dir="/mnt/galaxy/home/schavan/projects/bulk_rna_seq/RUN_Meloe1_celllines/"
sample_sheet = paste0(run_dir,"SampleSheet.txt") #Must have columns of sample, condition, batch, cancertype, provider

kal_dir = paste0(run_dir,"QUANT")
pca_raw_pdf = paste0(run_dir,"pca_raw.pdf")
pca_dds_pdf = paste0(run_dir,"pca_dds.pdf")
pca_limma_pdf = paste0(run_dir,"pca_adj.pdf")
gene_level_counts_adj  = paste0(run_dir,"gene_level_log2normcounts_adj.tsv")
gene_level_counts_adj_fil  = paste0(run_dir,"gene_level_log2normcounts_adj_fil.tsv")
gene_level_counts_tximport  = paste0(run_dir,"gene_level_rawcounts_kallisto_1.tsv")
gene_level_abundance_tximport = paste0(run_dir,"gene_level_abundanceTPM_kallisto_2.tsv")
kallisto_tpm_genelevel = paste0(run_dir,"gene_level_TPM_kallisto_2_with_genenames.tsv")
gene_level_scaled_tximport = paste0(run_dir,"gene_level_lengthScaledTPMCounts_kallisto_3.tsv")
deseq2_ncounts = paste0(run_dir,"deseq2_normalizedCnt_nobatch.csv")
boxplots_deseq2_pdf = paste0(run_dir,"deseq2_normalizedCnt.boxplots.nobatch.pdf")

### Edit above block to appropriate paths/filenames
###-----------------------------------------------------------------------------------------

kal_dir
sample_sheet

#s2c <- read.table(sample_sheet, header = FALSE, stringsAsFactors=FALSE)
s2c=fread(sample_sheet)
#s2c <- dplyr::select(s2c, sample = V1, condition = V2, batch = V3 , cancertype =V5)
s2c <- dplyr::mutate(s2c, path = paste0(kal_dir,"/",sample,".",condition))
s2c

files <- file.path(s2c$path, "abundance.h5")
names(files) <- paste0(s2c$sample,".",s2c$condition)
names(files)
all(file.exists(files))

for(i in 1:length(files))
{
    if(!file.exists(files[i])){
    print(paste0("File missing!-->",files[i]))}
}

tx2gene = fread("~/projects/bulk_rna_seq/HPVDetection/HPV-Hybrid/Hybrid_Kallisto_CustomMeloe1/transcripts_to_genes_w_hpv_genes_meloe1.txt",header=FALSE)
tx2gene = mutate(tx2gene, TXNAME=V1, GENEID=V2, SYMBOL=V3) %>% select(TXNAME, GENEID)
head(tx2gene)
print(dim(tx2gene))
#TXNAME            GENEID


###Store data in files
txi.kallisto.scale <- tximport(files, type = "kallisto", tx2gene = tx2gene, countsFromAbundance="lengthScaledTPM")
write.table(txi.kallisto.scale$counts, gene_level_scaled_tximport, row.names = F, sep="\t",quote=F,append=F)
dim(txi.kallisto.scale$counts)

txi.kallisto.tpm <- tximport(files, type = "kallisto", tx2gene = tx2gene) # countsFromAbundance = "no")
write.table(txi.kallisto.tpm$abundance, gene_level_abundance_tximport, row.names = F, sep="\t",quote=F,append=F)

txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene)
write.table(txi.kallisto$counts, gene_level_counts_tximport, row.names = F, sep="\t",quote=F,append=F)

### Add in gene names from the original tx2gene
tx2gene = fread("~/projects/bulk_rna_seq/HPVDetection/HPV-Hybrid/Hybrid_Kallisto_CustomMeloe1/transcripts_to_genes_w_hpv_genes_meloe1.txt",header=FALSE)
kallisto_tpm_with_genenames = txi.kallisto.tpm$abundance
rownames(kallisto_tpm_with_genenames) <- tx2gene$V3[match(rownames(kallisto_tpm_with_genenames), tx2gene$V2)]
head(kallisto_tpm_with_genenames)

### Make it a dataframe
df_kallisto_tpm_with_genenames <- cbind(Gene=rownames(kallisto_tpm_with_genenames),
                                     data.frame(kallisto_tpm_with_genenames, row.names=NULL))
head(df_kallisto_tpm_with_genenames)
df_kallisto_tpm_with_genenames$ENST = plyr::mapvalues(df_kallisto_tpm_with_genenames$Gene, tx2gene$V3, tx2gene$V1)
head(df_kallisto_tpm_with_genenames)
write.csv(df_kallisto_tpm_with_genenames, file=kallisto_tpm_genelevel, row.names = F);

names(txi.kallisto)
names(txi.kallisto$infReps)
df_txi_kal = as.data.frame(txi.kallisto$counts)
# transpose
  t_df_txi_kal <- transpose(df_txi_kal)
# get row and colnames in order
  colnames(t_df_txi_kal) <- rownames(df_txi_kal)
  rownames(t_df_txi_kal) <- colnames(df_txi_kal)

res.pca <- prcomp(t_df_txi_kal, scale = FALSE)
pca_scree = fviz_eig(res.pca)
pca_ind = fviz_pca_ind(res.pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)
pdf(pca_raw_pdf); pca_scree;pca_ind; #pca_var#pca_bi
dev.off()



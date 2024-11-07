library(readxl)
setwd("#####/project_directory")  # Set to project working directory
interactome_data <- read_xlsx("#####/interactome_map.xlsx",
                              sheet = 'Oligo interactome', skip = 2)[,1:6]  

enhancers <- read_xlsx("#####/interactome_map.xlsx",
                       sheet = 'Oligo enhancers', skip = 2,
                       col_names = c('chr', 'start', 'end'))

# Load gene expression data
load("#####/gene_expression_data.Rdata")
gene_expression <- as.data.frame(pbmc[["RNA"]]@counts)
gene_expression <- gene_expression[rowSums(gene_expression != 0) > 100,]

# Normalize and scale gene expression using Seurat
library(Seurat)
gexpr_dropout <- CreateSeuratObject(counts = gene_expression)
gexpr_dropout <- NormalizeData(gexpr_dropout)
gexpr_dropout <- ScaleData(gexpr_dropout, features = rownames(gexpr_dropout))

# Impute gene expression to address dropout
library(Rmagic)
sMatrix <- Matrix(data=as.matrix(gexpr_dropout@assays$RNA@data), sparse=TRUE)
gexpr_imputed <- magic(t(as.matrix(gexpr_dropout@assays$RNA@data)))

# Select cell type-specific gene expression (Oligodendrocytes)
table(pbmc@meta.data$broad.cell.type)
gexpr <- gexpr_imputed$result[rownames(pbmc@meta.data[pbmc@meta.data$broad.cell.type == "Oli",]),]
gexpr <- t(gexpr[,log10(colSums(gexpr)+1) > 1])

# Step 1: Identify chromatin interactions
df1 <- scGRNom_interaction(interactome_data, enhancers)
head(df1[sample(nrow(df1), 20),])

# Step 2: Infer transcription factor binding sites
df2 <- scGRNom_getTF(df1)
head(df2, 1)

# Step 3: Predict TF-target genes with elastic net regression
chromatin_access_regions <- read_xlsx("#####/chromatin_regions.xlsx", sheet = 'Feature Binarization Peaks', skip = 16)
open_chrom_regions <- chromatin_access_regions[chromatin_access_regions$Oligodendrocytes == 1,]
open_chrom_regions <- data.frame(na.omit(open_chrom_regions[,c('hg38_Chromosome', 'hg38_Start', 'hg38_Stop')]))
head(open_chrom_regions)

df3 <- scGRNom_getNt(df = df2, gexpr = gexpr, open_chrom = open_chrom_regions, extension_bps = 2000)
head(df3)

save(df2, df3, file = "#####/output_directory/Oli_tfbsnet.Rdata")

# Load existing TFBS network data
load("#####/output_directory/TFBS_network_data.Rdata")
ad_gwas <- data.frame(read.table("#####/gwas_data.txt"))[,c('CHR', 'BP', 'SNP')]
head(ad_gwas)

disease_genes <- c('PELP1', 'BIN1', 'MED11', 'ARRB2', 'PLCG2', 'PSMB6', 'PSMC3', 'SC5D', 'SORL1', 'C3', 'DENND1C', 
                   'CMIP', 'SUPT3H', 'RUNX2', 'SLC39A13', 'BCL7C', 'USP6NL', 'TMEM170B')

disease_genes <- scGRNom_disGenes(df = df3, gwas_snps = ad_gwas, extension_bps = 1000)
head(disease_genes)

# Analyze coefficient values in predicted interactions
df3ucoef <- df3[!duplicated(df3$coef),]
df3ucoef_top <- df3ucoef[abs(df3ucoef$coef) > 3,]
write.csv(df3ucoef_top, file = "#####/output_directory/top_tfbsnet_results.csv", row.names = FALSE)

# Statistical analysis on gene overlaps
q <- length(intersect(df3ucoef_top$TG, MT$X3))
m <- length(intersect(df2$gene, MT$X3))
phyper(q-1, m, length(unique(df2$gene)) - m, length(unique(df3ucoef_top$TG)), lower.tail = FALSE)

# GO and KEGG enrichment analysis
library(AnnotationHub)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(GOplot)

ego_gene <- bitr(unique(df3ucoef_MT$TG), fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = org.Hs.eg.db)

ego <- enrichGO(OrgDb="org.Hs.eg.db", keyType="ENSEMBL", gene=ego_gene$ENSEMBL, ont="BP", 
                pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.2, readable=TRUE)
dotplot(ego, showCategory=10, title="Top 10 GO Enrichment")
barplot(ego, showCategory=20, title="GO Enrichment")

# KEGG enrichment analysis
kegg_gene <- bitr(ego_gene$ENSEMBL, fromType="ENSEMBL", toType="UNIPROT", OrgDb=org.Hs.eg.db)
ekk <- enrichKEGG(gene=kegg_gene$UNIPROT, keyType="uniprot", organism='hsa', pvalueCutoff=0.05, 
                  pAdjustMethod="BH", qvalueCutoff=0.2)
dotplot(ekk, font.size=8, title='KEGG Enrichment')
test <- data.frame(ekk)
cnetplot(ekk, categorySize="pvalue", foldChange=ego_gene$ENSEMBL, colorEdge=TRUE, cex_label_gene=0.8)
heatplot(ekk)

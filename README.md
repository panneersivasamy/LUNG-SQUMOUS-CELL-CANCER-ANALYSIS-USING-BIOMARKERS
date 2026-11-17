# LUNG-SQUMOUS-CELL-CANCER-ANALYSIS-USING-BIOMARKERS
# =====================================================================
# Prognostic biomarkers in Squamous Cell Lung Carcinoma (LUSC) - Phase I
# Data   : TCGA-LUSC, STAR - Counts (via TCGAbiolinks)
# Method : DESeq2 (Tumor vs Normal), padj < 0.05 and |log2FC| > 1
# Run section by section (Ctrl + Enter). Do NOT source() the whole file
# unless the .rds files already exist (section 1 then skips the download).
# =====================================================================

# ---------------------------------------------------------------------
# 0. Packages
# ---------------------------------------------------------------------
# One-time install (uncomment if needed):
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("TCGAbiolinks", "DESeq2", "EnhancedVolcano"),
#                      ask = FALSE, update = FALSE)
# install.packages("pheatmap")

library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)

# ---------------------------------------------------------------------
# 1. Download and prepare TCGA-LUSC (skipped if already saved)
# ---------------------------------------------------------------------
if (file.exists("TCGA_LUSC_Data.rds")) {
  data <- readRDS("TCGA_LUSC_Data.rds")
} else {
  options(timeout = 600)
  httr::set_config(httr::config(connecttimeout = 120, timeout = 600, http_version = 0))

  query <- GDCquery(
    project       = "TCGA-LUSC",
    data.category = "Transcriptome Profiling",
    data.type     = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )

  # Retry loop for unstable networks (already downloaded files are skipped)
  for (i in 1:30) {
    dl <- try(GDCdownload(query, method = "api", files.per.chunk = 3), silent = TRUE)
    if (!inherits(dl, "try-error")) break
    message("Download retry ", i, " ...")
    Sys.sleep(15)
  }

  data <- GDCprepare(query)
  saveRDS(data, "TCGA_LUSC_Data.rds")
}

class(data)                 # RangedSummarizedExperiment
dim(data)                   # 60660 x 562
table(data$sample_type)     # 511 Primary Tumor, 51 Solid Tissue Normal

# ---------------------------------------------------------------------
# 2. Preprocessing and sample_info
# ---------------------------------------------------------------------
counts <- assay(data, "unstranded")          # raw counts

sample_info <- data.frame(
  row.names = colnames(counts),
  condition = ifelse(data$sample_type == "Solid Tissue Normal", "Normal", "Tumor")
)
sample_info$condition <- factor(sample_info$condition, levels = c("Normal", "Tumor"))
table(sample_info$condition)                 # Normal 51, Tumor 511

counts <- counts[rowSums(counts) >= 10, ]    # remove low-expression genes
dim(counts)                                  # 55299 x 562

# ---------------------------------------------------------------------
# 3. Differential expression (DESeq2)  - takes 10 to 30 minutes
# ---------------------------------------------------------------------
if (file.exists("dds.rds")) {
  dds <- readRDS("dds.rds")
} else {
  dds <- DESeqDataSetFromMatrix(countData = round(counts),
                                colData   = sample_info,
                                design    = ~ condition)
  dds <- DESeq(dds)
  saveRDS(dds, "dds.rds")
}

# NOTE: keep the name "res" only for the DESeq2 result
res <- results(dds, contrast = c("condition", "Tumor", "Normal"))
results_df <- as.data.frame(res)
write.csv(results_df, "DEG_Results.csv")
saveRDS(results_df, "results_df.rds")

significant_genes <- subset(results_df, padj < 0.05 & abs(log2FoldChange) > 1)
upregulated       <- subset(significant_genes, log2FoldChange > 1)
downregulated     <- subset(significant_genes, log2FoldChange < -1)
nrow(upregulated)       # 12350
nrow(downregulated)     # 5336

write.csv(upregulated,       "Upregulated_Genes.csv")
write.csv(downregulated,     "Downregulated_Genes.csv")
write.csv(significant_genes, "Significant_Genes.csv")

# ---------------------------------------------------------------------
# 4. Volcano plot
# ---------------------------------------------------------------------
gene_names <- rowData(data)$gene_name[match(rownames(results_df), rownames(data))]

volcano <- function() {
  EnhancedVolcano(results_df, lab = gene_names,
                  x = "log2FoldChange", y = "padj",
                  pCutoff = 0.05, FCcutoff = 1,
                  title = "LUSC: Tumor vs Normal")
}

print(volcano())                                   # shows in Plots tab
png("volcano_plot.png", width = 1200, height = 900, res = 150)
print(volcano())
dev.off()

# ---------------------------------------------------------------------
# 5. Heatmap (top 50 genes by adjusted p-value)
# ---------------------------------------------------------------------
vsd <- vst(dds, blind = FALSE)                     # 3 to 10 minutes

top50 <- head(order(results_df$padj), 50)
annotation_df <- data.frame(row.names = colnames(vsd),
                            condition = colData(dds)$condition)

pheatmap(assay(vsd)[top50, ], annotation_col = annotation_df,
         show_colnames = FALSE, show_rownames = FALSE, scale = "row")

png("heatmap_top50.png", width = 1400, height = 1000, res = 150)
pheatmap(assay(vsd)[top50, ], annotation_col = annotation_df,
         show_colnames = FALSE, show_rownames = FALSE, scale = "row")
dev.off()

# =====================================================================
# NEXT STEPS (not run yet) - GO/KEGG enrichment and survival analysis
# =====================================================================

# ---------------------------------------------------------------------
# 6. GO and KEGG enrichment
# ---------------------------------------------------------------------
# BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"),
#                      ask = FALSE, update = FALSE)
# library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)
#
# sig_ids  <- sub("\\..*", "", rownames(significant_genes))
# gene_ids <- bitr(sig_ids, fromType = "ENSEMBL", toType = "ENTREZID",
#                  OrgDb = org.Hs.eg.db)
#
# go_result <- enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db,
#                       ont = "BP", pAdjustMethod = "BH",
#                       pvalueCutoff = 0.05, readable = TRUE)
# png("GO_dotplot.png", width = 1400, height = 1000, res = 150)
# print(dotplot(go_result, showCategory = 15)); dev.off()
#
# kegg_result <- enrichKEGG(gene = gene_ids$ENTREZID, organism = "hsa")  # needs internet
# png("KEGG_dotplot.png", width = 1400, height = 1000, res = 150)
# print(dotplot(kegg_result, showCategory = 15)); dev.off()

# ---------------------------------------------------------------------
# 7. Survival analysis (Cox + Kaplan-Meier)
# ---------------------------------------------------------------------
# install.packages("survminer")
# library(survival); library(survminer)
#
# clinical <- as.data.frame(colData(data))
# clinical <- clinical[clinical$sample_type == "Primary Tumor", ]
# clinical <- clinical[!duplicated(clinical$patient), ]
# clinical$time   <- ifelse(clinical$vital_status == "Dead",
#                           clinical$days_to_death, clinical$days_to_last_follow_up)
# clinical$status <- ifelse(clinical$vital_status == "Dead", 1, 0)
# clinical <- clinical[!is.na(clinical$time) & clinical$time > 0, ]
#
# vsd_mat  <- assay(vsd)
# gene_map <- rowData(data)$gene_name
# names(gene_map) <- rownames(data)
#
# top <- head(rownames(significant_genes)[order(significant_genes$padj)], 200)
# cox_list <- lapply(top, function(ens) {
#   df <- clinical
#   df$expr <- vsd_mat[ens, rownames(df)]
#   cx <- summary(coxph(Surv(time, status) ~ expr, data = df))
#   data.frame(ensembl = ens, gene = gene_map[ens],
#              HR = cx$coefficients[1, "exp(coef)"],
#              p  = cx$coefficients[1, "Pr(>|z|)"])
# })
# cox_df <- do.call(rbind, cox_list)
# cox_df$padj <- p.adjust(cox_df$p, method = "BH")
# cox_df <- cox_df[order(cox_df$p), ]
# write.csv(cox_df, "Cox_results.csv", row.names = FALSE)
#
# plot_km <- function(sym) {
#   ens <- names(gene_map)[gene_map == sym][1]
#   df <- clinical
#   df$expr  <- vsd_mat[ens, rownames(df)]
#   df$group <- ifelse(df$expr > median(df$expr), "High", "Low")
#   fit <- survfit(Surv(time, status) ~ group, data = df)
#   ggsurvplot(fit, data = df, pval = TRUE, risk.table = TRUE, title = sym)
# }
# png("KM_plot.png", width = 1400, height = 1000, res = 150)
# print(plot_km(cox_df$gene[1])); dev.off()
file.exists("volcano_plot.png")
file.exists("heatmap_top50.png")
shell.exec("heatmap_top50.png")
graphics.off()
png("heatmap_top50.png", width = 1400, height = 1000, res = 150)
pheatmap(assay(vsd)[top50, ], annotation_col = annotation_df,
         show_colnames = FALSE, show_rownames = FALSE, scale = "row")
dev.off()
shell.exec("heatmap_top50.png")

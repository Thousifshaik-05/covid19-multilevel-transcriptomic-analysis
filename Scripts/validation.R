setwd("C://B.TECH//covid19-gene-expression-analysis//validation")

# ═══════════════════════════════════════════════════════════════════════
# VALIDATION OF DESEQ2 FINDINGS WITH EDGER AND LIMMA-VOOM
# Dataset: GSE152418 (COVID-19)
# ═══════════════════════════════════════════════════════════════════════

# ── Install packages if needed ─────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")
BiocManager::install("limma")
BiocManager::install("DESeq2")
BiocManager::install("ggvenn")

# ── Load libraries ─────────────────────────────────────────────────────
library(edgeR)
library(limma)
library(DESeq2)
library(ggplot2)
library(ggvenn)

# ── Set paths ──────────────────────────────────────────────────────────
results_path <- "C:/B.TECH/covid19-gene-expression-analysis/validation/val_res/"
dir.create(results_path, showWarnings = FALSE)

# ── Load raw counts ────────────────────────────────────────────────────
raw_counts <- read.csv(paste0(data_path, "GSE152418_p20047_Study1_RawCounts.txt"),
                       sep = "\t", row.names = 1)

cat("Raw counts dimensions:", dim(raw_counts), "\n")
head(raw_counts[1:3, 1:5])

# ── Load trait data ────────────────────────────────────────────────────
trait_data <- read.csv("C:/B.TECH/covid19-gene-expression-analysis/WGCNA/trait_data.csv")
head(trait_data)

# ── Prepare sample groups ──────────────────────────────────────────────
# Match sample order
raw_counts <- raw_counts[, trait_data$Sample]
group      <- factor(trait_data$Condition,
                     levels = c(0, 1),
                     labels = c("Control", "COVID"))
cat("Group table:\n")
print(table(group))

# ── Filter low count genes ─────────────────────────────────────────────
keep       <- rowSums(raw_counts >= 10) >= 3
raw_counts <- raw_counts[keep, ]
cat("Genes after filtering:", nrow(raw_counts), "\n")


#EDGER analysis

# ═══════════════════════════════════════════════════════════════
# EDGER
# ═══════════════════════════════════════════════════════════════
cat("\nRunning edgeR...\n")

dge        <- DGEList(counts = raw_counts, group = group)
dge        <- calcNormFactors(dge)
design     <- model.matrix(~group)
dge        <- estimateDisp(dge, design)
fit_edger  <- glmQLFit(dge, design)
qlf        <- glmQLFTest(fit_edger, coef = 2)

edger_results <- topTags(qlf, n = Inf)$table
edger_results$gene <- rownames(edger_results)

# Significant DEGs
edger_sig <- edger_results[edger_results$FDR < 0.05 &
                             abs(edger_results$logFC) >= 1, ]
cat("edgeR significant DEGs:", nrow(edger_sig), "\n")

# Save
write.csv(edger_results,
          paste0(results_path, "edgeR_All_Results.csv"),
          row.names = FALSE)
write.csv(edger_sig,
          paste0(results_path, "edgeR_Significant_DEGs.csv"),
          row.names = FALSE)
cat("✓ edgeR results saved\n")



#LIMMA-VOOM analysis

# ═══════════════════════════════════════════════════════════════
# LIMMA-VOOM
# ═══════════════════════════════════════════════════════════════
cat("\nRunning limma-voom...\n")

dge_limma  <- DGEList(counts = raw_counts, group = group)
dge_limma  <- calcNormFactors(dge_limma)
design     <- model.matrix(~group)
v          <- voom(dge_limma, design, plot = FALSE)
fit_limma  <- lmFit(v, design)
fit_limma  <- eBayes(fit_limma)

limma_results <- topTable(fit_limma, coef = 2, n = Inf)
limma_results$gene <- rownames(limma_results)

# Significant DEGs
limma_sig <- limma_results[limma_results$adj.P.Val < 0.05 &
                             abs(limma_results$logFC) >= 1, ]
cat("limma-voom significant DEGs:", nrow(limma_sig), "\n")

# Save
write.csv(limma_results,
          paste0(results_path, "limma_All_Results.csv"),
          row.names = FALSE)
write.csv(limma_sig,
          paste0(results_path, "limma_Significant_DEGs.csv"),
          row.names = FALSE)
cat("✓ limma-voom results saved\n")


#LOAD DESEQ2 AND COMPARE

# ═══════════════════════════════════════════════════════════════
# CONCORDANCE ANALYSIS
# ═══════════════════════════════════════════════════════════════

# Load your existing DESeq2 results
deseq2_all <- read.csv("C:/B.TECH/covid19-gene-expression-analysis/Results/All_DEG_results.csv")
deseq2_sig <- deseq2_all[!is.na(deseq2_all$padj) &
                           deseq2_all$padj < 0.05 &
                           abs(deseq2_all$log2FoldChange) >= 1, ]
cat("DESeq2 significant DEGs:", nrow(deseq2_sig), "\n")

# Gene lists from each method
deseq2_genes <- deseq2_sig$gene
edger_genes  <- edger_sig$gene
limma_genes  <- limma_sig$gene

cat("DESeq2 genes:", length(deseq2_genes), "\n")
cat("edgeR genes:",  length(edger_genes),  "\n")
cat("limma genes:",  length(limma_genes),  "\n")

# ── Concordance ────────────────────────────────────────────────
all3    <- Reduce(intersect, list(deseq2_genes, edger_genes, limma_genes))
d_e     <- intersect(deseq2_genes, edger_genes)
d_l     <- intersect(deseq2_genes, limma_genes)
e_l     <- intersect(edger_genes,  limma_genes)

cat("\n── Concordance Report ──\n")
cat("DESeq2 ∩ edgeR:          ", length(d_e),  "\n")
cat("DESeq2 ∩ limma:          ", length(d_l),  "\n")
cat("edgeR  ∩ limma:          ", length(e_l),  "\n")
cat("All 3 methods agree:     ", length(all3), "\n")
cat("Concordance rate (all3): ",
    round(length(all3) / length(deseq2_genes) * 100, 1), "%\n")

# Save concordant genes
concordant_df <- data.frame(gene = all3)
write.csv(concordant_df,
          paste0(results_path, "Concordant_DEGs_All3Methods.csv"),
          row.names = FALSE)
cat("✓ Concordant genes saved\n")


#VENN DIAGRAMS

# ── Venn diagram ───────────────────────────────────────────────
install.packages("ggvenn")
library(ggvenn)

venn_plot <- ggvenn(venn_list,
                    fill_color    = c("#e74c3c", "#3498db", "#2ecc71"),
                    stroke_size   = 0.5,
                    set_name_size = 5,
                    text_size     = 4) +
  labs(title = "Concordance of DEGs: DESeq2 vs edgeR vs limma-voom") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"))

ggsave(paste0(results_path, "Venn_DEG_Concordance.png"),
       plot   = venn_plot,
       width  = 8,
       height = 6,
       dpi    = 300,
       bg     = "white")
cat("✓ Venn diagram saved\n")

#CONCORDANCE PLOT

# Flip edgeR and limma logFC to match DESeq2 direction
edger_results$logFC <- -edger_results$logFC
limma_results$logFC <- -limma_results$logFC

# ── Log2FC correlation plot ────────────────────────────────────
# Merge DESeq2 and edgeR on common genes
common_de <- merge(deseq2_all[, c("gene", "log2FoldChange")],
                   edger_results[, c("gene", "logFC")],
                   by = "gene")
common_de <- merge(common_de,
                   limma_results[, c("gene", "logFC")],
                   by = "gene")
colnames(common_de) <- c("gene", "DESeq2", "edgeR", "limma")

# DESeq2 vs edgeR
scatter1 <- ggplot(common_de, aes(x = DESeq2, y = edgeR)) +
  geom_point(alpha = 0.3, size = 0.8, color = "#3498db") +
  geom_smooth(method = "lm", color = "red") +
  theme_bw() +
  labs(title = "Log2FC Correlation: DESeq2 vs edgeR",
       x = "DESeq2 log2FC", y = "edgeR logFC")

# DESeq2 vs limma
scatter2 <- ggplot(common_de, aes(x = DESeq2, y = limma)) +
  geom_point(alpha = 0.3, size = 0.8, color = "#2ecc71") +
  geom_smooth(method = "lm", color = "red") +
  theme_bw() +
  labs(title = "Log2FC Correlation: DESeq2 vs limma-voom",
       x = "DESeq2 log2FC", y = "limma logFC")

ggsave(paste0(results_path, "Scatter_DESeq2_vs_edgeR.png"),
       plot = scatter1, width = 6, height = 5, dpi = 300)
ggsave(paste0(results_path, "Scatter_DESeq2_vs_limma.png"),
       plot = scatter2, width = 6, height = 5, dpi = 300)
cat("✓ Scatter plots saved\n")

# Correlation coefficients
cat("\nCorrelation coefficients:\n")
cat("DESeq2 vs edgeR:", round(cor(common_de$DESeq2, common_de$edgeR), 3), "\n")
cat("DESeq2 vs limma:", round(cor(common_de$DESeq2, common_de$limma), 3), "\n")
cat("edgeR  vs limma:", round(cor(common_de$edgeR,  common_de$limma), 3), "\n")

# ── Libraries ──────────────────────────────────────────────────────────
library(GEOquery)
library(ggplot2)
library(tidyr)
library(dplyr)
library(org.Hs.eg.db)

# ── Results folder ─────────────────────────────────────────────────────
results_path <- "C:/B.TECH/covid19-gene-expression-analysis/Data/Results/"
dir.create(results_path, showWarnings = FALSE)

# ── Load metadata ──────────────────────────────────────────────────────
gse      <- getGEO("GSE152418", GSEMatrix = TRUE)
metadata <- pData(gse[[1]])

clinical <- metadata[, c("geo_accession",
                         "disease state:ch1",
                         "severity:ch1",
                         "gender:ch1",
                         "days_post_symptom_onset:ch1")]
colnames(clinical) <- c("sample", "disease_state", "severity", "gender", "days_post_onset")
clinical$sample_name <- gsub("-", ".", metadata$title)

cat("Disease state:\n"); print(table(clinical$disease_state))
cat("Severity:\n");      print(table(clinical$severity))

# ── Hub genes ──────────────────────────────────────────────────────────
hub_genes_final <- c("CDC45", "SPAG5", "TTK", "CCNA2", "BIRC5",
                     "BUB1B", "PBK", "CDK1", "KIF2C", "MELK")

hub_ensembl <- mapIds(org.Hs.eg.db,
                      keys    = hub_genes_final,
                      column  = "ENSEMBL",
                      keytype = "SYMBOL")
hub_ensembl <- na.omit(hub_ensembl)
print(hub_ensembl)

# ── Load normalized counts ─────────────────────────────────────────────
norm_counts <- read.csv("C:/B.TECH/covid19-gene-expression-analysis/WGCNA/normalized_counts.csv")
rownames(norm_counts) <- norm_counts$X
norm_counts <- norm_counts[, -which(colnames(norm_counts) == "X")]
cat("norm_counts dimensions:", dim(norm_counts), "\n")

# ── Extract hub gene expression ────────────────────────────────────────
hub_expr <- norm_counts[rownames(norm_counts) %in% hub_ensembl, ]
cat("hub_expr dimensions:", dim(hub_expr), "\n")  # should be 10 x 34

hub_expr_t <- as.data.frame(t(hub_expr))
ensembl_to_symbol <- setNames(names(hub_ensembl), hub_ensembl)
colnames(hub_expr_t) <- ensembl_to_symbol[colnames(hub_expr_t)]
hub_expr_t$sample_name <- rownames(hub_expr_t)
cat("hub_expr_t columns:", colnames(hub_expr_t), "\n")

# ── Merge ──────────────────────────────────────────────────────────────
combined <- merge(hub_expr_t, clinical, by = "sample_name")
cat("Combined dimensions:", dim(combined), "\n")
cat("Combined columns:", colnames(combined), "\n")

# Save
write.csv(combined,
          paste0(results_path, "HubGenes_Clinical_Combined.csv"),
          row.names = FALSE)
cat("✓ Combined data saved\n")

# ── Severity boxplot ───────────────────────────────────────────────────
hub_cols    <- c("sample_name", "severity", names(hub_ensembl))
hub_subset  <- combined[, hub_cols]
hub_long    <- pivot_longer(hub_subset,
                            cols      = all_of(names(hub_ensembl)),
                            names_to  = "gene",
                            values_to = "expression")
hub_long$severity <- factor(hub_long$severity,
                            levels = c("Healthy","Moderate","Severe","ICU","Convalescent"))

severity_plot <- ggplot(hub_long, aes(x = severity, y = expression, fill = severity)) +
  geom_boxplot() +
  facet_wrap(~gene, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Hub Gene Expression Across Severity Groups",
       x = "Severity", y = "Normalized Expression") +
  scale_fill_manual(values = c("Healthy"      = "#2ecc71",
                               "Moderate"     = "#f1c40f",
                               "Severe"       = "#e67e22",
                               "ICU"          = "#e74c3c",
                               "Convalescent" = "#3498db"))
ggsave(paste0(results_path, "Severity_HubGenes.png"),
       plot = severity_plot, width = 14, height = 10, dpi = 300)
write.csv(hub_long,
          paste0(results_path, "Severity_HubGenes_Data.csv"),
          row.names = FALSE)
cat("✓ Severity plot saved\n")

# ── Disease state boxplot ──────────────────────────────────────────────
disease_cols   <- c("sample_name", "disease_state", names(hub_ensembl))
disease_subset <- combined[, disease_cols]
disease_long   <- pivot_longer(disease_subset,
                               cols      = all_of(names(hub_ensembl)),
                               names_to  = "gene",
                               values_to = "expression")
disease_long$disease_state <- factor(disease_long$disease_state,
                                     levels = c("Healthy","COVID-19","Convalescent"))

disease_plot <- ggplot(disease_long, aes(x = disease_state, y = expression, fill = disease_state)) +
  geom_boxplot() +
  facet_wrap(~gene, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Hub Gene Expression: COVID-19 vs Healthy vs Convalescent",
       x = "Disease State", y = "Normalized Expression") +
  scale_fill_manual(values = c("Healthy"      = "#2ecc71",
                               "COVID-19"     = "#e74c3c",
                               "Convalescent" = "#3498db"))
ggsave(paste0(results_path, "DiseaseState_HubGenes.png"),
       plot = disease_plot, width = 14, height = 10, dpi = 300)
write.csv(disease_long,
          paste0(results_path, "DiseaseState_HubGenes_Data.csv"),
          row.names = FALSE)
cat("✓ Disease state plot saved\n")

# ── Severity comparison (Moderate vs Severe vs ICU) ───────────────────
severity_only <- combined[combined$severity %in% c("Moderate","Severe","ICU"), ]
severity_only$severity <- factor(severity_only$severity,
                                 levels = c("Moderate","Severe","ICU"))
cat("Samples per severity group:\n")
print(table(severity_only$severity))

sev_subset  <- severity_only[, hub_cols]
severity_long <- pivot_longer(sev_subset,
                              cols      = all_of(names(hub_ensembl)),
                              names_to  = "gene",
                              values_to = "expression")

# Kruskal-Wallis test
kruskal_results <- data.frame(gene = names(hub_ensembl), p_value = NA)
for(i in seq_along(names(hub_ensembl))){
  gene      <- names(hub_ensembl)[i]
  gene_data <- severity_long[severity_long$gene == gene, ]
  test      <- kruskal.test(expression ~ severity, data = gene_data)
  kruskal_results$p_value[i] <- round(test$p.value, 4)
}
kruskal_results$significance <- ifelse(kruskal_results$p_value < 0.05,
                                       "Significant", "Not significant")
write.csv(kruskal_results,
          paste0(results_path, "Severity_Kruskal_Statistics.csv"),
          row.names = FALSE)
cat("\nKruskal-Wallis results:\n")
print(kruskal_results)

severity_compare_plot <- ggplot(severity_long,
                                aes(x = severity, y = expression, fill = severity)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 2) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
  facet_wrap(~gene, scales = "free_y", ncol = 5) +
  theme_bw() +
  theme(axis.text.x   = element_text(angle = 45, hjust = 1),
        strip.text     = element_text(face = "bold"),
        legend.position = "bottom") +
  labs(title = "Hub Gene Expression: Moderate vs Severe vs ICU",
       x = "Severity", y = "Normalized Expression") +
  scale_fill_manual(values = c("Moderate" = "#f1c40f",
                               "Severe"   = "#e67e22",
                               "ICU"      = "#e74c3c"))
ggsave(paste0(results_path, "Severity_Comparison_Plot.png"),
       plot = severity_compare_plot, width = 16, height = 8, dpi = 300)
cat("✓ Severity comparison plot saved\n")

# ── Cross-disease comparison (COVID vs Flu vs RSV) ────────────────────
gse_flu  <- getGEO("GSE111368", GSEMatrix = TRUE)
flu_expr <- exprs(gse_flu[[1]])
flu_meta <- pData(gse_flu[[1]])

gse_rsv  <- getGEO("GSE34205", GSEMatrix = TRUE)
rsv_expr <- exprs(gse_rsv[[1]])
rsv_meta <- pData(gse_rsv[[1]])

# Flu hub gene extraction
flu_features <- fData(gse_flu[[1]])
flu_probes   <- flu_features[flu_features$Symbol %in% hub_genes_final,
                             c("ID", "Symbol")]
cat("Flu probes found:\n"); print(flu_probes)

flu_expr_hub <- as.data.frame(t(flu_expr[rownames(flu_expr) %in% flu_probes$ID, ]))
colnames(flu_expr_hub) <- flu_probes$Symbol[match(colnames(flu_expr_hub), flu_probes$ID)]
flu_expr_hub$sample_name <- rownames(flu_expr_hub)
flu_expr_hub$disease     <- ifelse(flu_meta$`flu_type:ch1` == "HC", "Healthy_Flu", "Influenza")

# RSV hub gene extraction
rsv_features <- fData(gse_rsv[[1]])
cat("RSV feature columns:\n"); print(colnames(rsv_features))



# ── RSV hub gene extraction ────────────────────────────────────────────

colnames(rsv_features)
head(rsv_features[, 1:5])
# ── RSV hub gene extraction ────────────────────────────────────────────
rsv_probes <- rsv_features[rsv_features$`Gene Symbol` %in% hub_genes_final,
                           c("ID", "Gene Symbol")]
cat("RSV probes found:\n")
print(rsv_probes)

rsv_expr_hub <- as.data.frame(t(rsv_expr[rownames(rsv_expr) %in% rsv_probes$ID, ]))
colnames(rsv_expr_hub) <- rsv_probes$`Gene Symbol`[match(colnames(rsv_expr_hub), 
                                                         rsv_probes$ID)]
rsv_expr_hub$sample_name <- rownames(rsv_expr_hub)
rsv_expr_hub$disease     <- ifelse(rsv_meta$`infection:ch1` == "none", "Healthy_RSV",
                                   ifelse(rsv_meta$`infection:ch1` == "rsv",  "RSV", "Influenza_RSV"))

cat("RSV samples:", nrow(rsv_expr_hub), "\n")
print(head(rsv_expr_hub[, c("sample_name", "disease")]))




# ── Prepare COVID data ─────────────────────────────────────────────────
covid_expr_hub <- combined[, c("sample_name", names(hub_ensembl))]
covid_expr_hub$disease <- "COVID-19"

# ── Normalize column names across all datasets ─────────────────────────
# Keep only common hub genes across all 3 datasets
common_genes <- Reduce(intersect, list(
  names(hub_ensembl),
  colnames(flu_expr_hub),
  colnames(rsv_expr_hub)
))
cat("Common hub genes across all datasets:", length(common_genes), "\n")
print(common_genes)

# ── Subset each dataset to common genes ───────────────────────────────
covid_sub <- covid_expr_hub[, c(common_genes, "disease")]
flu_sub   <- flu_expr_hub[,  c(common_genes, "disease")]
rsv_sub   <- rsv_expr_hub[,  c(common_genes, "disease")]

# ── Combine all three datasets ─────────────────────────────────────────
all_diseases <- rbind(covid_sub, flu_sub, rsv_sub)
cat("Combined dataset dimensions:", dim(all_diseases), "\n")
table(all_diseases$disease)

# ── Reshape to long format ─────────────────────────────────────────────
all_long <- pivot_longer(all_diseases,
                         cols      = all_of(common_genes),
                         names_to  = "gene",
                         values_to = "expression")

# Set disease order
all_long$disease <- factor(all_long$disease,
                           levels = c("Healthy_RSV", "Healthy_Flu",
                                      "COVID-19", "Influenza",
                                      "Influenza_RSV", "RSV"))

# ── Cross-disease comparison plot ──────────────────────────────────────
cross_disease_plot <- ggplot(all_long,
                             aes(x = disease, y = expression, fill = disease)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1) +
  facet_wrap(~gene, scales = "free_y", ncol = 5) +
  theme_bw() +
  theme(axis.text.x    = element_text(angle = 45, hjust = 1),
        strip.text      = element_text(face = "bold"),
        legend.position = "bottom") +
  labs(title = "Hub Gene Expression: COVID-19 vs Influenza vs RSV",
       x     = "Disease",
       y     = "Expression") +
  scale_fill_manual(values = c("Healthy_RSV"   = "#2ecc71",
                               "Healthy_Flu"   = "#27ae60",
                               "COVID-19"      = "#e74c3c",
                               "Influenza"     = "#e67e22",
                               "Influenza_RSV" = "#f39c12",
                               "RSV"           = "#3498db"))

# ── Save plot ──────────────────────────────────────────────────────────
ggsave(paste0(results_path, "CrossDisease_Comparison_Plot.png"),
       plot   = cross_disease_plot,
       width  = 18,
       height = 10,
       dpi    = 300)
cat("✓ Cross-disease plot saved\n")

# ── Save data ──────────────────────────────────────────────────────────
write.csv(all_long,
          paste0(results_path, "CrossDisease_HubGenes_Data.csv"),
          row.names = FALSE)
cat("✓ Cross-disease data saved\n")

# ── Statistical test across diseases ──────────────────────────────────
cat("\nKruskal-Wallis test per gene across diseases:\n")
kruskal_cross <- data.frame(gene = common_genes, p_value = NA)
for(i in seq_along(common_genes)){
  gene      <- common_genes[i]
  gene_data <- all_long[all_long$gene == gene, ]
  test      <- kruskal.test(expression ~ disease, data = gene_data)
  kruskal_cross$p_value[i] <- round(test$p.value, 4)
}
kruskal_cross$significance <- ifelse(kruskal_cross$p_value < 0.05,
                                     "Significant", "Not significant")
print(kruskal_cross)

write.csv(kruskal_cross,
          paste0(results_path, "CrossDisease_Kruskal_Statistics.csv"),
          row.names = FALSE)
cat("✓ Cross-disease statistics saved\n")





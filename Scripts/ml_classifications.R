setwd("C://B.TECH//covid19-gene-expression-analysis//ML_Classsification")


# ═══════════════════════════════════════════════════════════════════════
# MACHINE LEARNING CLASSIFICATION — COVID-19 vs HEALTHY
# Using Hub Gene Signature (10 genes)
# Methods: Random Forest + SVM
# ═══════════════════════════════════════════════════════════════════════

# ── Install packages ───────────────────────────────────────────────────
install.packages("randomForest")
install.packages("e1071")          # SVM
install.packages("caret")
install.packages("pROC")
install.packages("ggplot2")

# ── Load libraries ─────────────────────────────────────────────────────
library(randomForest)
library(e1071)
library(caret)
library(pROC)
library(ggplot2)
library(org.Hs.eg.db)

# ── Paths ──────────────────────────────────────────────────────────────
results_path <- "C:/B.TECH/covid19-gene-expression-analysis/ML_Classsification/"
dir.create(results_path, showWarnings = FALSE)

# ── Load normalized counts ─────────────────────────────────────────────
norm_counts <- read.csv("C:/B.TECH/covid19-gene-expression-analysis/WGCNA/normalized_counts.csv")
rownames(norm_counts) <- norm_counts$X
norm_counts <- norm_counts[, -which(colnames(norm_counts) == "X")]

# ── Hub genes ──────────────────────────────────────────────────────────
hub_genes_final <- c("CDC45", "SPAG5", "TTK", "CCNA2", "BIRC5",
                     "BUB1B", "PBK", "CDK1", "KIF2C", "MELK")

hub_ensembl <- mapIds(org.Hs.eg.db,
                      keys    = hub_genes_final,
                      column  = "ENSEMBL",
                      keytype = "SYMBOL")
hub_ensembl <- na.omit(hub_ensembl)

# ── Extract hub gene expression ────────────────────────────────────────
hub_expr <- norm_counts[rownames(norm_counts) %in% hub_ensembl, ]
hub_expr_t <- as.data.frame(t(hub_expr))
ensembl_to_symbol <- setNames(names(hub_ensembl), hub_ensembl)
colnames(hub_expr_t) <- ensembl_to_symbol[colnames(hub_expr_t)]

cat("Expression matrix:", dim(hub_expr_t), "\n")

# ── Load labels ────────────────────────────────────────────────────────
trait_data <- read.csv("C:/B.TECH/covid19-gene-expression-analysis/WGCNA/trait_data.csv")
trait_data$Sample <- gsub("-", ".", trait_data$Sample)

# Match sample order
hub_expr_t$Sample <- rownames(hub_expr_t)
ml_data <- merge(hub_expr_t, trait_data, by = "Sample")
ml_data$Sample <- NULL

# Convert condition to factor
ml_data$Condition <- factor(ml_data$Condition,
                            levels = c(0, 1),
                            labels = c("Healthy", "COVID"))
cat("Class distribution:\n")
print(table(ml_data$Condition))
cat("Final ML dataset dimensions:", dim(ml_data), "\n")

# ── Train/Test Split ───────────────────────────────────────────────────
set.seed(42)
train_idx   <- createDataPartition(ml_data$Condition, p = 0.7, list = FALSE)
train_data  <- ml_data[train_idx, ]
test_data   <- ml_data[-train_idx, ]

cat("Training samples:", nrow(train_data), "\n")
cat("Testing samples:",  nrow(test_data),  "\n")

# ── Cross-validation setup ─────────────────────────────────────────────
ctrl <- trainControl(
  method          = "cv",
  number          = 5,
  classProbs      = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = TRUE
)

# ═══════════════════════════════════════════════════════════════════════
# RANDOM FOREST
# ═══════════════════════════════════════════════════════════════════════
cat("\nTraining Random Forest...\n")


# Reload caret
library(caret)

# Now run again
set.seed(42)
rf_model <- caret::train(
  Condition ~ .,
  data       = train_data,
  method     = "rf",
  metric     = "ROC",
  trControl  = ctrl,
  importance = TRUE
)
# Predictions
rf_pred       <- predict(rf_model, test_data)
rf_pred_prob  <- predict(rf_model, test_data, type = "prob")

# Confusion matrix
cat("\nRandom Forest Results:\n")
rf_cm <- confusionMatrix(rf_pred, test_data$Condition)
print(rf_cm)

# ROC
rf_roc <- roc(test_data$Condition,
              rf_pred_prob$COVID,
              levels = c("Healthy", "COVID"))
cat("Random Forest AUC:", round(auc(rf_roc), 3), "\n")

# Variable importance
# Add overall importance as mean of both classes
rf_importance$Overall <- rowMeans(rf_importance[, c("Healthy", "COVID")])

# Sort
rf_importance <- rf_importance[order(rf_importance$Overall, decreasing = TRUE), ]
cat("\nFeature Importance:\n")
print(rf_importance)

# ── Feature Importance Plot ────────────────────────────────────────────
imp_plot <- ggplot(rf_importance,
                   aes(x = reorder(gene, Overall), y = Overall)) +
  geom_bar(stat = "identity", fill = "#e74c3c") +
  coord_flip() +
  theme_bw() +
  labs(title = "Random Forest — Feature Importance",
       x     = "Gene",
       y     = "Importance Score")

ggsave(paste0(results_path, "RF_Feature_Importance.png"),
       plot   = imp_plot,
       width  = 8,
       height = 6,
       dpi    = 300)
cat("✓ Feature importance plot saved\n")
# ═══════════════════════════════════════════════════════════════════════
# SVM
# ═══════════════════════════════════════════════════════════════════════
cat("\nTraining SVM...\n")

library(caret)
set.seed(42)

svm_model <- caret::train(
  Condition ~ .,
  data      = train_data,
  method    = "svmRadial",
  metric    = "ROC",
  trControl = ctrl,
  preProcess = c("center", "scale")
)

# Predictions
svm_pred      <- predict(svm_model, test_data)
svm_pred_prob <- predict(svm_model, test_data, type = "prob")

# Confusion matrix
cat("\nSVM Results:\n")
svm_cm <- confusionMatrix(svm_pred, test_data$Condition)
print(svm_cm)

# ROC
svm_roc <- roc(test_data$Condition,
               svm_pred_prob$COVID,
               levels = c("Healthy", "COVID"))
cat("SVM AUC:", round(auc(svm_roc), 3), "\n")

# ═══════════════════════════════════════════════════════════════════════
# PLOTS
# ═══════════════════════════════════════════════════════════════════════

# ── ROC Curve ─────────────────────────────────────────────────────────
png(paste0(results_path, "ROC_Curve_RF_SVM.png"),
    width = 800, height = 600)

plot(rf_roc,  col = "#e74c3c", lwd = 2,
     main = "ROC Curves: Random Forest vs SVM")
plot(svm_roc, col = "#3498db", lwd = 2, add = TRUE)
legend("bottomright",
       legend = c(paste0("Random Forest (AUC = ", round(auc(rf_roc),  3), ")"),
                  paste0("SVM           (AUC = ", round(auc(svm_roc), 3), ")")),
       col    = c("#e74c3c", "#3498db"),
       lwd    = 2)
dev.off()
cat("✓ ROC curve saved\n")

# ── Feature Importance Plot ────────────────────────────────────────────
imp_plot <- ggplot(rf_importance,
                   aes(x = reorder(gene, Overall), y = Overall)) +
  geom_bar(stat = "identity", fill = "#e74c3c") +
  coord_flip() +
  theme_bw() +
  labs(title = "Random Forest — Feature Importance",
       x     = "Gene",
       y     = "Importance Score")

ggsave(paste0(results_path, "RF_Feature_Importance.png"),
       plot   = imp_plot,
       width  = 8,
       height = 6,
       dpi    = 300)
cat("✓ Feature importance plot saved\n")

# ── Concordance summary table ──────────────────────────────────────────
results_summary <- data.frame(
  Method   = c("Random Forest", "SVM"),
  AUC      = c(round(auc(rf_roc),  3),
               round(auc(svm_roc), 3)),
  Accuracy = c(round(rf_cm$overall["Accuracy"],  3),
               round(svm_cm$overall["Accuracy"], 3)),
  Sensitivity = c(round(rf_cm$byClass["Sensitivity"],  3),
                  round(svm_cm$byClass["Sensitivity"], 3)),
  Specificity = c(round(rf_cm$byClass["Specificity"],  3),
                  round(svm_cm$byClass["Specificity"], 3))
)

print(results_summary)
write.csv(results_summary,
          paste0(results_path, "ML_Results_Summary.csv"),
          row.names = FALSE)
cat("✓ ML results summary saved\n")





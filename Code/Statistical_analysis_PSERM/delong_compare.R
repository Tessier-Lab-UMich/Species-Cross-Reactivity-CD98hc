#!/usr/bin/env Rscript
# Perform DeLong's test to compare the two ROC curves of two metrics
# -----------------------------
# Load required package
# -----------------------------
if(!require(pROC)) {
  install.packages("pROC")
  library(pROC)
}

# -----------------------------
# Read input CSV file
# -----------------------------
data <- read.csv("./variants_PSERM.csv")
#data <- read.csv("./variants_PSERM_ER_D1.csv")


labels <- data$labels
scores1 <- data$FRE_D1
#scores1 <- data$PSSM_D1
#scores1 <- data$ER_D1
scores2 <- data$PSERM_D1

# -----------------------------
# Compute ROC curves for the two scoring metrics
# -----------------------------
roc1 <- roc(labels, scores1)
roc2 <- roc(labels, scores2)

auc1 <- auc(roc1)
auc2 <- auc(roc2)

cat("AUC 1 =", auc1, "\n")
cat("95% CI 1 =", ci.auc(roc1)[1], "-", ci.auc(roc1)[3], "\n\n")

cat("AUC 2 =", auc2, "\n")
cat("95% CI 2 =", ci.auc(roc2)[1], "-", ci.auc(roc2)[3], "\n\n")

# -----------------------------
# Perform DeLong's test to compare the two ROC curves
# -----------------------------
delong_test <- roc.test(roc1, roc2, method="delong")
cat("DeLong Z =", delong_test$statistic, "\n")
cat("p-value =", delong_test$p.value, "\n\n")






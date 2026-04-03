## Compute AUC and plot ROC curve for predicting PLN winner
## Predictors: log(N/D) and MAC (mean absolute correlation)
## Outcome: winner == "PLN" (binary, 1 = PLN wins)

library(pROC)

data_path <- "/projects/genomic-ml/da2343/PLN/pln_eval/out/count_prediction_results_retained20_with_stats.csv"
out_dir   <- "/projects/genomic-ml/da2343/PLN/pln_eval/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

df <- read.csv(data_path, stringsAsFactors = FALSE)
df$pln_wins  <- as.integer(df$winner == "PLN")
df$log_nd    <- log(df$N_over_D)

# ── Logistic regression models ────────────────────────────────────────────────

fit_nd      <- glm(pln_wins ~ log_nd,        data = df, family = binomial)
fit_mac     <- glm(pln_wins ~ mac,            data = df, family = binomial)
fit_both    <- glm(pln_wins ~ log_nd + mac,   data = df, family = binomial)

pred_nd     <- predict(fit_nd,   type = "response")
pred_mac    <- predict(fit_mac,  type = "response")
pred_both   <- predict(fit_both, type = "response")

roc_nd      <- roc(df$pln_wins, pred_nd,   quiet = TRUE)
roc_mac     <- roc(df$pln_wins, pred_mac,  quiet = TRUE)
roc_both    <- roc(df$pln_wins, pred_both, quiet = TRUE)

auc_nd   <- round(as.numeric(auc(roc_nd)),   3)
auc_mac  <- round(as.numeric(auc(roc_mac)),  3)
auc_both <- round(as.numeric(auc(roc_both)), 3)

cat("AUC (N/D only):       ", auc_nd,   "\n")
cat("AUC (MAC only):       ", auc_mac,  "\n")
cat("AUC (N/D + MAC):      ", auc_both, "\n")

# ── Plot ROC curves ───────────────────────────────────────────────────────────

out_png <- file.path(out_dir, "roc_auc_nd_mac.png")
png(out_png, width = 700, height = 620, res = 110)

plot(roc_nd,   col = "#2166AC", lwd = 2.5, main = "ROC curves: predicting PLN win",
     xlab = "False Positive Rate (1 - Specificity)",
     ylab = "True Positive Rate (Sensitivity)")
plot(roc_mac,  col = "#D6604D", lwd = 2.5, add = TRUE)
plot(roc_both, col = "#1A9641", lwd = 2.5, add = TRUE)
abline(a = 0, b = 1, lty = 2, col = "grey60")

legend("bottomright", bty = "n", lwd = 2.5,
       col  = c("#2166AC", "#D6604D", "#1A9641"),
       legend = c(
         sprintf("N/D only  (AUC = %.3f)", auc_nd),
         sprintf("MAC only  (AUC = %.3f)", auc_mac),
         sprintf("N/D + MAC (AUC = %.3f)", auc_both)
       ))

dev.off()
cat("Saved ROC plot to:", out_png, "\n")

# ── Model summaries ───────────────────────────────────────────────────────────
cat("\n--- N/D only ---\n");   print(summary(fit_nd))
cat("\n--- MAC only ---\n");   print(summary(fit_mac))
cat("\n--- N/D + MAC ---\n");  print(summary(fit_both))

results_df <-read.csv(csvname)
pdf(pdfname)
#Bin size equal to alpha
bin_width <- 0.05
breaks_seq <- seq(from = floor(min(results_df$p_value, na.rm = TRUE) / bin_width) * bin_width,
                  to = ceiling(max(results_df$p_value, na.rm = TRUE) / bin_width) * bin_width,
                  by = bin_width)
                
#Plotting
hist(results_df$p_value,
     prob = TRUE,
     main = NULL,
     xlab = "p-value",
     ylab = "Density",
     col = "#C7A9ABff",
     border = "black",
     breaks = breaks_seq,
     ylim = c(0, 3))
abline(h = 1, col = "#280003ff", lwd = 3, lty = "dashed")
dev.off()

#FPR
sig_p <- results_df$p_value[results_df$p_value < 0.025 | results_df$p_value > 0.975]
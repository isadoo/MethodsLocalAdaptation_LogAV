ds <- read.csv(csvname)

pdf(pdfname)


#Bins with size of alpha!
bin_width <- 0.05
breaks_seq <- seq(from = floor(min(ds$S_Value, na.rm = TRUE) / bin_width) * bin_width,
                  to = ceiling(max(ds$S_Value, na.rm = TRUE) / bin_width) * bin_width,
                  by = bin_width)
#Plotting
hist(ds$S_Value,
     prob = TRUE,
     main = NULL,
     xlab = "S-value",
     ylab = "Density",
     col = "#7B5154",
     border = "black",
     breaks = 20,
     xlim = c(0, 1))
     
abline(v = 0.5, col = "#280003ff", lwd = 3, lty = "dashed")

#Printing
dev.off()

#FPR
sig_S <- ds$S_Value[ds$S_Value < 0.2 | ds$S_Value > 0.8]
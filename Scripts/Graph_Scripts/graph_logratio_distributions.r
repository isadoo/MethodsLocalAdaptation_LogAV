# Load required packages
library(ggplot2)
library(dplyr)
library(readr)
library(purrr)
library(tidyr)

csv_directory <- "Chapter1/Simulations_Chapter1"  

#List all CSV files in the directory
csv_files <- list.files(path = csv_directory, pattern = "\\.csv$", full.names = TRUE)

#to extract log_ratio and clean label
extract_log_ratio <- function(file) {
  data <- read.csv(file)
  #Extract label from filename (everything after 'Lava')
  label <- sub(".*Lava", "Lava", basename(file))
  data.frame(log_ratio = data$log_ratio, file = label)
}

#Combine all into one dataframe
combined_data <- do.call(rbind, lapply(csv_files, extract_log_ratio))

#unique labels
labels <- unique(combined_data$file)


combined_data$file <- factor(combined_data$file, levels = labels)


plotmy <- ggplot(combined_data, aes(x = file, y = `log_ratio`, fill = file)) +
  geom_violin(trim = FALSE) +
  theme_minimal(base_size = 14) +
  labs(title = "Distribution of log_ratio",
       x = "Population and environment",
       y = "log ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

pdf("Chapter1/Simulations_Chapter1/graph_logratio_distributions.pdf", width = 10, height = 6)
print(plotmy)
dev.off()


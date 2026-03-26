library(hierfstat)
library(JGTeach)
library(Matrix)
library(gaston)
library(MCMCglmm)
library(brms)
library(ggplot2)
library(parallel)

library(ggplot2)
library(dplyr)
library(readr)
library(purrr)
library(tidyr)

simulation_directory <- "Chapter1/Simulations_Chapter1/quantinemo_" 




#Default breeding parameters
ns <- 5  # Number of sires
nd <- 5  # Number of dams
no <- 2  # Number of offspring per pair
n_loci_qtl <- 100  # Number of QTLs
maplength_quanti <- 50  # Map length for drop along ped


process_replicate <- function(file_path, pop_structure) {
  if (!file.exists(file_path)) {
    return(NULL)
  }
  
  
  sim_quanti <- read.fstat(fname = file_path)
  dos_quanti <- biall2dos(sim_quanti[, -1])
  
 
  traits <- rowSums((dos_quanti - 1) * 0.2) + rnorm(nrow(dos_quanti))
  
  
  if (grepl("IM", pop_structure)) {
    n_pops <- 8
    ind_per_pop <- 10
  } else {  # SS
    n_pops <- 20
    ind_per_pop <- 10
  }
  
  pop <- rep(1:n_pops, each = ind_per_pop)
  
  return(data.frame(traits = traits, pop = pop))
}


process_simulation_type <- function(sim_type) {
  cat("Processing", sim_type, "...\n")
  
  
  if (grepl("IM", sim_type)) {
    generation <- 500
    n_pops <- 8
  } else {  # SS
    generation <- 5000
    n_pops <- 20
  }
  
  
  if (sim_type == "IM_Neutral") {
    
    base_dir <- paste0("Chapter1/Simulations_Chapter1/quantinemo_", sim_type, "/")
    all_traits <- list()
    
    for (rep in 1:500) {
      rep_str <- sprintf("%03d", rep)
      file_path <- paste0(base_dir, "quanti_trait_g", generation, "_r", rep_str, ".dat")
      
      rep_data <- process_replicate(file_path, sim_type)
      if (!is.null(rep_data)) {
        rep_data$replicate <- rep
        all_traits[[rep]] <- rep_data
      }
    }
  } else {
    
    base_dir <- paste0("Chapter1/Simulations_Chapter1/quantinemo_", sim_type, "/")
    all_traits <- list()
    
    for (rep in 1:500) {
      rep_dir <- paste0(base_dir, "quantinemo_", sim_type, "_rep", rep, "/")
      file_path <- paste0(rep_dir, "quanti_trait_g", generation, ".dat")
      
      rep_data <- process_replicate(file_path, sim_type)
      if (!is.null(rep_data)) {
        rep_data$replicate <- rep
        all_traits[[rep]] <- rep_data
      }
    }
  }
  
  
  if (length(all_traits) > 0) {
    combined_data <- do.call(rbind, all_traits)
    combined_data$sim_type <- sim_type
    return(combined_data)
  } else {
    return(NULL)
  }
}

#Process all simulation types
sim_types <- c("IM_Neutral", "IM_omega10", "IM_omega22", "IM_omega50",
               "SS_Neutral", "SS_20pop_omega10", "SS_20pop_omega22", "SS_20pop_omega50")

all_sim_results <- list()


for (sim_type in sim_types) {
  result <- process_simulation_type(sim_type)
  all_sim_results[[sim_type]] <- result
}


create_boxplot <- function(data, sim_type) {
  
  pop_stats <- data %>%
    group_by(pop) %>%
    summarise(
      mean_trait = mean(traits, na.rm = TRUE),
      median_trait = median(traits, na.rm = TRUE),
      q25 = quantile(traits, 0.25, na.rm = TRUE),
      q75 = quantile(traits, 0.75, na.rm = TRUE),
      ci_lower = quantile(traits, 0.025, na.rm = TRUE),
      ci_upper = quantile(traits, 0.975, na.rm = TRUE),
      .groups = 'drop'
    )
  
  
  par(mar = c(5, 4, 4, 2) + 0.1)
  
 
  trait_list <- split(data$traits, data$pop)
  
  
  bp <- boxplot(trait_list, 
                main = paste("Trait Distribution -", sim_type),
                xlab = "Population",
                ylab = "Trait Value",
                col = "lightblue",
                border = "darkblue",
                las = 1)
  
  
  for (i in 1:nrow(pop_stats)) {
    lines(c(i, i), c(pop_stats$ci_lower[i], pop_stats$ci_upper[i]), 
          col = "red", lwd = 2)
    lines(c(i-0.1, i+0.1), c(pop_stats$ci_lower[i], pop_stats$ci_lower[i]), 
          col = "red", lwd = 2)
    lines(c(i-0.1, i+0.1), c(pop_stats$ci_upper[i], pop_stats$ci_upper[i]), 
          col = "red", lwd = 2)
  }
  
  
  legend("topright", 
         legend = c("Boxplot", "95% CI"), 
         col = c("darkblue", "red"), 
         lwd = c(1, 2),
         bty = "n")
  
  return(pop_stats)
}

pdf("Chapter1/Simulations_Chapter1/trait_distribution_plots.pdf", width = 12, height = 8)

par(mfrow = c(2, 4)) 

summary_stats <- list()

for (sim_type in names(all_sim_results)) {
  if (!is.null(all_sim_results[[sim_type]])) {
    cat("Creating plot for", sim_type, "...\n")
    stats <- create_boxplot(all_sim_results[[sim_type]], sim_type)
    summary_stats[[sim_type]] <- stats
  }
}
dev.off()

par(mfrow = c(1, 1))


cat("\nSummary Statistics:\n")
for (sim_type in names(summary_stats)) {
  cat("\n", sim_type, ":\n")
  print(summary_stats[[sim_type]])
}


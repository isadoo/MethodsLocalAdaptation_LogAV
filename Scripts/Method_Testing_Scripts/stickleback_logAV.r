library(dplyr)
library(tidyr)
library(hierfstat)
library(brms)
library(Matrix)
source("Chapter1/applying_logAV/functions_applying_logav.r")

#1. get the data
trait_data <- read.csv("Chapter1/applying_logAV/stickleback_karhunen/karhunen2014_phenotypes_standardized.csv")
genetic_P <- read.csv("Chapter1/applying_logAV/stickleback_karhunen/karhunen2014_genotypes_dosage.csv")
pedigree <- read.csv("Chapter1/applying_logAV/stickleback_karhunen/karhunen2014_pedigree_standardized.csv")

#2. transform our genetic data into the right format
#Create dosage matrix: 
dosage_P <- genetic_P %>% 
  select(-pop)  
dosage_P <- as.matrix(dosage_P[,-1]) 
population_P <- genetic_P$pop

head(pedigree)
head(trait_data)
str(genetic_P)
str(dosage_P)

population_id <- data.frame(id = pedigree$id, "population_id" = pedigree$dam_pop)

#3. Get thetaP and M
coancestries <- ThetaP_TheM(genetic_data_parents = dosage_P,
                            genotyped_parent_populations = population_P,
                            population_individual_id = population_id,
                            pedigree = pedigree)

#add proper names to the matrices
#For The.M: use individual IDs from pedigree
rownames(coancestries$The.M) <- pedigree$id
colnames(coancestries$The.M) <- pedigree$id

#For Theta.P: use unique population names
unique_populations <- unique(population_P)
rownames(coancestries$Theta.P) <- unique_populations
colnames(coancestries$Theta.P) <- unique_populations

cat("Added names to matrices:\n")
cat("The.M dimensions:", dim(coancestries$The.M), "with rownames:", head(rownames(coancestries$The.M)), "...\n")
cat("Theta.P dimensions:", dim(coancestries$Theta.P), "with rownames:", rownames(coancestries$Theta.P), "\n")

#4. we check if The.M and Theta.P are pd
eigT <- eigen(coancestries$Theta.P)
eigM <- eigen(coancestries$The.M)
min(eigT$values) 
min(eigM$values) 

#5. Now focus on trait data
#we check for NA in the trait_data
sum(is.na(trait_data))

#6.Create The.M per trait, such that we take out the individuals with NA from the The.M specific of that trait


unique(trait_data$trait_id)

#

#[1] "SL"          "body.depth"  "cp.length"   "feed"        "head.length"
#[6] "n.attack"    "pg.length"   "sex"         "t.fullout"   "t.orient"
#We will only consider the following traits:
#SL = "Standard Length"
#body.depth = "Body Depth"
#cp.length = "Caudal Peduncle Length"
#head.length = "Head Length"
#n.attack = "Aggression"
#pg.lenght = "Pelvic gridle Length"
#t.fullout = "Risk taking 1" = "time in seconds till the fish fully came out from a refuge"
#feed = "Rist taking 2" = "whether fish bit on food after a simulated attack within 300 sec"
#t.orient = "timet to orient in seconds"

#get list of all unique traits (excluding sex which we'll use as covariate)
traits_to_analyze <- c("SL", "body.depth", "cp.length", "head.length", "n.attack", "pg.length", "t.fullout", "feed", "t.orient")

#create trait-specific matrices and data
trait_matrices <- list()
trait_data_filtered_list <- list()

for (trait_name in traits_to_analyze) {
  cat("\n=== Processing trait:", trait_name, "===\n")
  
  #Filtering trait data to get individuals without NA for this specific trait
  trait_data_specific <- trait_data %>%
    filter(trait_id == trait_name, !is.na(trait))
  
  #get individual IDs that have data for this trait
  individuals_with_data <- trait_data_specific$individual
  
  #check which of these individuals are also in The.M matrix
  individuals_in_matrix <- individuals_with_data[individuals_with_data %in% rownames(coancestries$The.M)]
  
  cat("Total individuals with", trait_name, "data:", length(individuals_with_data), "\n")
  cat("Individuals also in The.M matrix:", length(individuals_in_matrix), "\n")
  
  if (length(individuals_in_matrix) < 5) {
    cat("Warning: Too few individuals for", trait_name, "- skipping\n")
    next
  }
  
  #IMPORTANT: Order individuals by their position in the original The.M matrix
  #This preserves the block diagonal structure by population
  original_order <- rownames(coancestries$The.M)
  individuals_in_matrix_ordered <- original_order[original_order %in% individuals_in_matrix]
  
  #create trait-specific The.M matrix by subsetting in the correct order
  The.M_trait <- coancestries$The.M[individuals_in_matrix_ordered, individuals_in_matrix_ordered]
  
  #check if matrix is PD
  eig_trait <- eigen(The.M_trait)
  cat("Minimum eigenvalue for", trait_name, ":", min(eig_trait$values), "\n")
  
  #make positive definite if needed
  if (min(eig_trait$values) < 0) {
    cat("Making matrix positive definite for", trait_name, "\n")
    pd_The.M_trait <- as.matrix(nearPD(The.M_trait)$mat)
    #restore the names after nearPD (which removes them)
    rownames(pd_The.M_trait) <- individuals_in_matrix_ordered
    colnames(pd_The.M_trait) <- individuals_in_matrix_ordered
  } else {
    pd_The.M_trait <- The.M_trait
  }
  
  #filter trait data to match the matrix individuals (in correct order)
  trait_data_filtered <- trait_data_specific %>%
    filter(individual %in% individuals_in_matrix_ordered) %>%
    #Order the trait data to match the matrix order
    arrange(match(individual, individuals_in_matrix_ordered))
  
  #check population structure using pedigree data
  trait_pops_from_pedigree <- trait_data_filtered %>%
    left_join(pedigree %>% select(id, dam_pop), by = c("individual" = "id"))
  
  if("dam_pop" %in% colnames(trait_pops_from_pedigree)) {
    pop_counts <- table(trait_pops_from_pedigree$dam_pop)
    cat("Population distribution for", trait_name, ":", paste(names(pop_counts), "=", pop_counts, collapse = ", "), "\n")
  } else {
    cat("Could not determine population distribution for", trait_name, "\n")
  }
  
  
  trait_matrices[[trait_name]] <- pd_The.M_trait
  trait_data_filtered_list[[trait_name]] <- trait_data_filtered
  
  cat("Successfully created matrix for", trait_name, "with", nrow(pd_The.M_trait), "individuals\n")
}

cat("\n=== Summary of trait matrices ===\n")
for (trait_name in names(trait_matrices)) {
  cat(trait_name, ": matrix size", nrow(trait_matrices[[trait_name]]), "x", ncol(trait_matrices[[trait_name]]), "\n")
}


#We will run logAV on these ones, and we can consider sex as a covariate

#run logAV analysis for each trait
logAV_results <- list()

for (trait_name in names(trait_matrices)) {
  cat("\n=== Running logAV for trait:", trait_name, "===\n")
  
  #Get the trait-specific matrix and data
  pd_The.M_trait <- trait_matrices[[trait_name]]
  trait_data_filtered <- trait_data_filtered_list[[trait_name]]
  
  #Add population information from pedigree
  trait_data_with_pop <- trait_data_filtered %>%
    left_join(pedigree %>% select(id, dam_pop), by = c("individual" = "id")) %>%
    rename(population = dam_pop)
  
  #Add sex information as covariate if available
  trait_data_with_sex <- trait_data_with_pop %>%
    left_join(
      trait_data %>% 
        filter(trait_id == "sex") %>% 
        select(individual, sex_value = trait),
      by = "individual"
    )
  
  #Set up formula covariates
  if (sum(!is.na(trait_data_with_sex$sex_value)) > 0) {
    formula_covariates <- "(1 | sex_value)"
    cat("Including sex as covariate for", trait_name, "\n")
  } else {
    formula_covariates <- NULL
    cat("No sex covariate available for", trait_name, "\n")
  }
  
  #run logAV analysis
  tryCatch({
    if (!is.null(formula_covariates)) {
      logAV_result <- logAV(Theta.P = coancestries$Theta.P, 
                           The.M = pd_The.M_trait, 
                           trait_dataframe = trait_data_with_sex,
                           column_population = "population",
                           column_individual = "individual", 
                           column_trait = "trait",
                           formula_covariates = formula_covariates)
    } else {
      logAV_result <- logAV(Theta.P = coancestries$Theta.P, 
                           The.M = pd_The.M_trait, 
                           trait_dataframe = trait_data_with_pop,
                           column_population = "population",
                           column_individual = "individual", 
                           column_trait = "trait")
    }
    
    logAV_results[[trait_name]] <- logAV_result
    cat("Successfully analyzed", trait_name, "\n")
    
  }, error = function(e) {
    cat("Error analyzing", trait_name, ":", e$message, "\n")
    logAV_results[[trait_name]] <- list(error = e$message)
  })
}


cat("\n=== SUMMARY OF RESULTS ===\n")
for (trait_name in names(logAV_results)) {
  if ("error" %in% names(logAV_results[[trait_name]])) {
    cat(trait_name, ": FAILED\n")
  } else {
    result <- logAV_results[[trait_name]]
    cat(sprintf("%-15s: Mean log ratio = %8.4f, P-value = %8.4f\n", 
                trait_name, result$BRMS_mean_log_ratio, result$p_value))
  }
}




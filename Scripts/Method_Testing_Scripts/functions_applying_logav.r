##### Functions ------------------------------------------------------------------------------------------------
ThetaP_TheM <- function(genetic_data_parents,
                                   genotyped_parent_populations, 
                                   genetic_data_F1 = NA, 
                                   population_individual_id = NA, 
                                   column_individual = "id", 
                                   column_population = "population_id", 
                                   pedigree = NA,
                                   all_parents_genotyped = FALSE) {



    parent_dosage <- genetic_data_parents
    

    # Load F1 data and get kinship ------------------------------------------------------ 
    if (!identical(pedigree, NA)) {
        kinship_F1 <- kinship_from_pedigree(pedigree) #internal function from package LAVA
        #Make sure individuals are arranged by their naming, which ideally follows the population naming
        #sorted_F1_names <- sort(rownames(kinship_F1))
        #kinship_F1 <- kinship_F1[sorted_F1_names, sorted_F1_names]

        pop_ids_F1 <- pedigree$dam_pop

        F1_id <- pedigree$id
    } else if (!identical(genetic_data_F1, NA)) {
        F1_dosage <- genetic_data_F1
       
        matching_matrix_F1 <- hierfstat::matching(F1_dosage) 
        kinship_F1 <- hierfstat::beta.dosage(matching_matrix_F1, MATCHING = TRUE)
        F1_id <- population_individual_id[,column_individual]
    } else {
        stop("Missing F1 data, either include pedigree or genetic data")
    }
   
    # ----------------------------------------------------------------------------

    # Determine parental population IDs and F1 population IDs (which includes P) --------------------
    if (!identical(pedigree, NA)) {
        population_names <- unique(pedigree$dam_pop)
        parent_pop_id <- genotyped_parent_populations #may differ from pedigree
        #Calculate population sizes for The.M
        dam_df <- pedigree %>%
            select(id = dam, population_id = dam_pop) %>%
            distinct()
        sire_df <- pedigree %>%
            select(id = sire, population_id = sire_pop) %>%
            distinct()
        id_df <- pedigree %>%
            select(id, population_id = dam_pop) %>%
            distinct()

        parent_pop_id_pedigree <- c(sire_df$population_id, dam_df$population_id)
        F1_pop_id <- id_df$population_id

        dam_per_pop <- dam_df %>%
                count(population_id)
        sire_per_pop <- sire_df %>% 
                count(population_id)
        F1_per_pop <- id_df %>%
                count(population_id)
        sire_per_pop <- sire_per_pop %>%
                rename(sire_n = n)
        dam_per_pop <- dam_per_pop %>%
            rename(dam_n = n)
        parent_per_pop <- full_join(sire_per_pop, dam_per_pop, by = "population_id") %>%
                mutate(
                    sire_n = replace_na(sire_n, 0),
                    dam_n = replace_na(dam_n, 0),
                    parent_n = sire_n + dam_n
                )
        
        population_sizes_F1 <- F1_per_pop$n
        population_sizes_P <- parent_per_pop$parent_n #On pedigree (may differ from genotyped)
        population_sizes <- population_sizes_P + population_sizes_F1 #On pedigree (may differ from genotyped)
    } else if (!identical(population_individual_id, NA)) {
        
            population_names <- unique(population_individual_id[, column_population])
            parent_pop_id <- genotyped_parent_populations
            F1_pop_id <- population_individual_id[, column_population] 
            #Calculate population sizes for The.M
            population_sizes_P <- as.data.frame(table(parent_pop_id))$Freq
            population_sizes_F1 <- as.data.frame(table(F1_pop_id))$Freq
            population_sizes <- population_sizes_P + population_sizes_F1
         
    } else {
            stop("Population data missing. \n
            Provide either population_individual_id \n")
    }
    

    cat("There are ", length(unique(parent_pop_id)), " populations. \n")

    if (exists("pop_ids_F1")) {
        if (!identical(pop_ids_F1,F1_pop_id)) {
            warning("Mismatch between ids in F1 genetic data or pedigree and population identification.\n")
        }
    } else {
        pop_ids_F1 <- F1_pop_id
    }



    #Theta_P calculation -----------------------------------------------------
    #Matching matrix and kinship for parents
    matching_matrix_parents <- hierfstat::matching(parent_dosage)
    kinship_parents <- hierfstat::beta.dosage(matching_matrix_parents, MATCHING = TRUE)
    fst_founders <- hierfstat::fs.dosage(matching_matrix_parents, pop = parent_pop_id, matching = TRUE)

    cat("Calculating Theta.P \n")
    min_Fst <- min(hierfstat::mat2vec(fst_founders$FsM))
    Theta_P <- (fst_founders$FsM - min_Fst) / (1 - min_Fst)
    cat("Theta.P calculated with dimensions", dim(Theta_P), "\n")
    # --------------------------------------------------------------------------


    # The.M calculation -------------------------------------------------------
    cat("Calculating The.M \n")
    
    # Check if population sizes are correct ----------------------
    total_individuals <- sum(population_sizes)
    total_individuals_F1 <- sum(population_sizes_F1)
    #cat("population sizes including F1 and Founders are ", population_sizes, "\n")
    #cat("total individuals are ", total_individuals, "\n")
    cat("population sizes of F1", population_sizes_F1, "\n")
    cat("total individuals are ", total_individuals_F1, "\n")

    

    #Calculate mean kinship per population of parental (founder) population
    unique_pops <- unique(F1_pop_id) #same number of populations in F1 and in P
    #Should be the same as population_names

    

    
    #Explaining following If statement:
    #If you don't have genetic info from the parents, 
    #then you cannot infer their relationship, 
    #so you have to assume it is zero 
    #when you use the pedigree to estimate relatedness in F1
    if(all_parents_genotyped == FALSE || !identical(pedigree, NA)){
        adjusted_kinship_F1 <- kinship_F1
    } else {
        #Adjust kinship for F1 individuals - we use the mean kinship of the founders to standardize M
        adjusted_kinship_F1 <- matrix(0,nrow=nrow(kinship_F1),ncol=nrow(kinship_F1))
        last_individual_counted <- 0
        mean_kinship_per_population <- sapply(unique_pops, function(pop) {
        pop_indices_P <- which(parent_pop_id == pop)
        mean(hierfstat::mat2vec(kinship_parents[pop_indices_P, pop_indices_P]))
        })
        
        for (pop in 0:(length(unique_pops) - 1)) {  #Version before May 14th 2025 allowed for dosage of phenotyped individuals to include founders. 
        #Find indices for this population's F1
            cat("pop: ", pop, "\n")

            #Adjusting F1 individuals' kinship matrix
            start_idx <- last_individual_counted + 1
            end_idx <- last_individual_counted + population_sizes_F1[pop+1]
            cat("start_idx: ", start_idx, "end_idx: ", end_idx, "\n")
                
            adjusted_kinship_F1[start_idx:end_idx, start_idx:end_idx] <- 
                (kinship_F1[start_idx:end_idx, start_idx:end_idx] - mean_kinship_per_population[pop+1]) / 
                (1 - mean_kinship_per_population[pop+1])

            #Update the last individual count tracker
            last_individual_counted <- end_idx
            }

    }

    
    The.M <- matrix(0, nrow = sum(population_sizes_F1), ncol = sum(population_sizes_F1))
    row.names(The.M) <- colnames(The.M) <- F1_id


    for (pop in unique_pops) {
    #Grab the indices for individuals from current population
        pop_indices <- which(pop_ids_F1 == pop)
        dim_current_population <- length(pop_indices)
        cat("There are", dim_current_population, "individuals in population", pop, "for calculating the M matrix\n")
        kin_block <- hierfstat::kinship2grm(adjusted_kinship_F1)[pop_indices, pop_indices]
        block_adjusted <- kin_block * (1 - Theta_P[pop, pop])
        The.M[pop_indices, pop_indices] <- block_adjusted
    }

    #The M must be positive definite
    eigenvalues <- eigen(The.M)$values
    if (any(eigenvalues < 0)) {
        cat("M matrix not positive definite. \n")
        cat("Minimum eigenvalue is ", min(eigenvalues), "\n")
    }

    cat("The.M calculated with dimensions ", dim(The.M), "\n")

    return(list(The.M = The.M, Theta.P = Theta_P))
}

logAV <- function(Theta.P, 
                The.M, 
                trait_dataframe,
                column_population = "population", 
                column_individual = "id", 
                column_trait = "trait",
                formula_covariates = NULL) {
  
  #Check input types and dimensions
  if (!is.matrix(Theta.P) || !is.matrix(The.M)) {
    stop("Theta.P and The.M must be matrices.")
  }
  
  if (!is.data.frame(trait_dataframe) || ncol(trait_dataframe) < 2) {
    stop("trait_dataframe must be a data frame with at least two columns (ID and trait values).")
  }
  
  #Identify populations per individual
  population_blocks_df <- counting_blocks_matrix(The.M)
  individuals_per_population_F1 <- population_blocks_df$rows
  number_of_blocks <- length(population_blocks_df$block) 
  pop_ids <- rep(1:number_of_blocks, individuals_per_population_F1[1:number_of_blocks]) 
  number_populations <- nrow(Theta.P)
  
  #if (number_of_blocks != number_populations) {
  #  warning(paste0("Mismatch between detected populations based on The.M matrix (", number_of_blocks,") and Theta.P dimensions (", number_populations, ")."))
  #}
  
  #Standardize trait data - FIX THE LIST ISSUE
  Y <- trait_dataframe[, column_trait]
  
  # Debug: Check what we're getting
  cat("Class of Y before processing:", class(Y), "\n")
  cat("Is Y a list?", is.list(Y), "\n")
  
  # Handle the list issue properly
  if (is.list(Y)) {
    Y <- unlist(Y)
  }
  
  # Ensure Y is numeric
  Y <- as.numeric(Y)
  
  # Remove any NAs
  valid_indices <- !is.na(Y)
  Y <- Y[valid_indices]
  
  # Filter the dataframe to match
  trait_dataframe <- trait_dataframe[valid_indices, ]
  
  # Debug: Check final Y
  cat("Class of Y after processing:", class(Y), "\n")
  cat("Length of Y:", length(Y), "\n")
  cat("Is Y a list after processing?", is.list(Y), "\n")
  
  # Standardize
  Y <- Y - mean(Y)
  var_Y <- var(Y)
  Y <- Y / sqrt(var_Y)
  
  #From VB = VA*2FST
  two.Theta.P <- 2 * Theta.P
  
  #Prepare data frame for modeling
  dat <- trait_dataframe
  dat$Y <- Y  # This should now be a numeric vector
  
  # Fix the pop and ind columns - they might also be lists
  pop_col <- trait_dataframe[,column_population]
  if (is.list(pop_col)) {
    pop_col <- unlist(pop_col)
  }
  dat$pop <- as.character(pop_col)
  
  ind_col <- trait_dataframe[,column_individual]
  if (is.list(ind_col)) {
    ind_col <- unlist(ind_col)
  }
  dat$ind <- as.character(ind_col)
  
  # Debug: Check all columns
  cat("Is dat$Y a list?", is.list(dat$Y), "\n")
  cat("Is dat$pop a list?", is.list(dat$pop), "\n")
  cat("Is dat$ind a list?", is.list(dat$ind), "\n")
  cat("Class of dat$pop:", class(dat$pop), "\n")
  cat("Class of dat$ind:", class(dat$ind), "\n")
  
  #Build the complete formula
  base_formula <- "Y ~ 1 + (1 | gr(pop, cov = two.Theta.P)) + (1 | gr(ind, cov = The.M))"
  
  if (!is.null(formula_covariates)) {
    formula_string <- paste0("Y ~ 1 + ", formula_covariates, " + (1 | gr(pop, cov = two.Theta.P)) + (1 | gr(ind, cov = The.M))")
  } else {
    formula_string <- base_formula
  }
  
  model_formula <- as.formula(formula_string)
  
  cat("Using formula:", formula_string, "\n")
  
  #Bayesian model with covariates
  brms_mf <- brm(model_formula, 
                 data = dat, 
                 data2 = list(two.Theta.P = two.Theta.P, The.M = The.M), 
                 family = gaussian(), 
                 chains = 8, 
                 cores = 4, 
                 iter = 3000, 
                 warmup = 1000, 
                 thin = 2)
  
  #variance components
  var_components <- lapply(VarCorr(brms_mf, summary = FALSE), function(x) x$sd^2)
  var_df <- as.data.frame(do.call(cbind, var_components))
  
  quant_med <- quantile(var_df$pop - var_df$ind, c(0.5, 0.025, 0.975))
  mean_diff <- mean(var_df$pop - var_df$ind)
  
  #Hypothesis testing
  hyp <- "sd_pop__Intercept^2 - sd_ind__Intercept^2 = 0"
  the_hyp <- hypothesis(brms_mf, hyp, class = NULL)
  
  #Posteriors: VA,B and VA,A
  post_samples <- as_draws_df(brms_mf, variable = c("sd_pop__Intercept", "sd_ind__Intercept"))
  post_samples$log_ratio <- log(post_samples$sd_pop__Intercept^2 / post_samples$sd_ind__Intercept^2)
  
  mean_log_ratio <- mean(post_samples$log_ratio)
  quant_log_ratio <- quantile(post_samples$log_ratio, probs = c(0.025, 0.975))
  
  p_value <- 2 * mean(sign(post_samples$log_ratio) != sign(median(post_samples$log_ratio)))
  
  results <- list(
    p_value = p_value,
    post_samples = post_samples,
    brms_model = brms_mf,  
    BRMS_test = the_hyp$hypothesis[2:5],
    BRMS_mean_diff = mean_diff,
    BRMS_test_median = quant_med["50%"],
    BRMS_test_median_lower = quant_med["2.5%"],
    BRMS_test_median_upper = quant_med["97.5%"],
    BRMS_mean_log_ratio = mean_log_ratio,
    BRMS_log_ratio_ci_lower = quant_log_ratio[1],
    BRMS_log_ratio_ci_upper = quant_log_ratio[2],
    formula_used = formula_string
  )
  
  return(results)
}

kinship_from_pedigree <- function(pedigree) {
  
  ids <- pedigree$id
  num_individuals <- length(ids)
  kinship_matrix <- matrix(0, nrow = num_individuals, ncol = num_individuals)
  rownames(kinship_matrix) <- ids
  colnames(kinship_matrix) <- ids
  
  #Self-relatedness = 1 
  diag(kinship_matrix) <- 1
  
  
  for (i in seq_len(num_individuals)) {
    for (j in seq_len(i - 1)) { 
      sire_i <- pedigree$sire[i]
      dam_i <- pedigree$dam[i]
      sire_j <- pedigree$sire[j]
      dam_j <- pedigree$dam[j]
      
      if (!is.na(sire_i) && !is.na(sire_j) && sire_i == sire_j) {
        kinship_matrix[i, j] <- kinship_matrix[j, i] <- kinship_matrix[i, j] + (1/4) #half sib
      }
      if (!is.na(dam_i) && !is.na(dam_j) && dam_i == dam_j) {
        kinship_matrix[i, j] <- kinship_matrix[j, i] <- kinship_matrix[i, j] + (1/4) #half sib
      }
      #If you're "halfsib" on both sides, then  you're full sib!
    }
  }
  
  return(kinship_matrix)
}

counting_blocks_matrix <- function(mat) {
  n <- nrow(mat)
  blocks <- list()
  visited <- rep(FALSE, n)  # Track visited rows
  block_count <- 0

  for (i in seq_len(n)) {
    if (!visited[i]) {
      #Find the start of a new block
      block_count <- block_count + 1
      block_rows <- which(mat[i, ] != 0)

      #Check for contiguous rows forming a block
      block_size <- 1
      while ((i + block_size) <= n && 
             all(mat[i + block_size, block_rows] != 0) &&
             all(mat[i + block_size, -block_rows] == 0)) {
        block_size <- block_size + 1
      }
      
      #Mark the rows as visited
      visited[i:(i + block_size - 1)] <- TRUE
      
      #Store block information
      blocks[[block_count]] <- list(block = block_count, rows = block_size, cols = length(block_rows))
    }
  }

  
  block_df <- do.call(rbind, lapply(blocks, as.data.frame))
  return(block_df)
}


######################################-------------------------------------------------------------------------------


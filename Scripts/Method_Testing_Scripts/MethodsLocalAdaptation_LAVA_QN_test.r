.libPaths("/work/FAC/FBM/DEE/jgoudet/default/R/lib")
library(hierfstat)
library(JGTeach)
library(Matrix)
library(gaston)
library(MCMCglmm)
library(brms)
library(ggplot2)
library(parallel)

#Getting data ########################################################################

args <- commandArgs(trailingOnly = TRUE)

#Arguments to the respective variables
replicate_number <- as.integer(args[1])
optima_typeCheck <- args[2]
if(optima_typeCheck == "Neutral"){optima_type = ""} else {optima_type = "Selection"}
generations <- as.integer(args[3])
Population_structure <- args[4]
Selective_or_Neutral <- args[5]
np <- as.integer(args[6])

#File paths and names
Simulation_directory <- paste0("/work/FAC/FBM/DEE/jgoudet/pop_fst/isa/Chapter1/QN_", 
                                Population_structure, "_", Selective_or_Neutral, optima_type, "/")
Neutral_file_name <- paste0("neutral_data_g", generations, "_r")
Quanti_file_name <- paste0("quanti_trait_g", generations, "_r")
file_csv_name <- paste0("MethodsLocalAdaptation_Lava_", Selective_or_Neutral, "_", Population_structure)

#Default breeding parameters
ns <- 5  # Number of sires
nd <- 5  # Number of dams
no <- 2  # Number of offspring per pair
n_replicates <- 100  # Number of replicates
n_loci_qtl <- 100  # Number of QTLs
maplength_quanti <- 50  # Map length for drop along ped
##################################################################################################

#Print out arguments either passed by user or default (available in output)
cat("optima_type:", optima_type, "\n")
cat("Population_structure:", Population_structure, "\n")
cat("Simulation_directory:", Simulation_directory, "\n")
cat("Neutral_file_name:", Neutral_file_name, "\n")
cat("Quanti_file_name:", Quanti_file_name, "\n")
cat("ns:", ns, "\n")
cat("nd:", nd, "\n")
cat("np:", np, "\n")
cat("no:", no, "\n")
cat("n_replicates:", n_replicates, "\n")


###############################################################################################

#Setting pedigree (the same for all replicates)####
#total number of founders
nft<-np*(ns+nd)
sire<-rep(1:ns,each=nd*no)
dam<-rep(ns+1:nd,each=no,ns)
sire<-rep(0:(np-1)*(nd+ns),each=(nd*ns*no))+sire
dam<-rep(0:(np-1)*(nd+ns),each=(nd*ns*no))+dam
# need to have the parents of founders as NA
sire<-c(rep(NA,nft),sire)
dam<-c(rep(NA,nft),dam)
nt<-length(sire)
nf<- -c(1:nft)
ni<-ns+nd
#create the pedigree, 
pedi<-data.frame(ind=1:nt,sire,dam)
###################################################


get_data <- function(dat, dat_quanti){
    dos<-biall2dos(dat[,-1])
    nl<-ncol(dos)

    #dos_quanti <- fstat2dos(dat_quanti[, -1])
    dos_quanti <- biall2dos(dat_quanti[, -1])
    nl_quanti <- ncol(dos_quanti)

    combined_nl <- nl + nl_quanti
    combined_founders <- array(0, dim = c(nft, combined_nl, 2))
    combined_dat <- cbind(dat[, -1], dat_quanti[, -1])

    tmp1<-as.matrix(combined_dat%/%10-1)
    tmp2<-as.matrix(combined_dat%%10-1)
    combined_founders[,,1]<-tmp1
    combined_founders[,,2]<-tmp2

    #combined_founders[, 1:nl, 1] <- as.matrix(dat[, -1] %/% 10 - 1)
    #combined_founders[, 1:nl, 2] <- as.matrix(dat[, -1] %% 10 - 1)
    #combined_founders[, (nl + 1):combined_nl, 1] <- as.matrix(dat_quanti[, -1] %/% 10 - 1)
    #combined_founders[, (nl + 1):combined_nl, 2] <- as.matrix(dat_quanti[, -1] %% 10 - 1)
    nfounders_combined <- rbind(combined_founders[, , 1], combined_founders[, , 2])

    genos.F1_combined <- drop.along.ped(pedi, founders.genotypes = nfounders_combined, nloc = combined_nl, maplength = 1000)

    #Extract the neutral and quantitative genotypes from the combined result
    genos.F1_neutral <- genos.F1_combined[, 1:nl, ]
    genos.F1_quanti <- genos.F1_combined[, (nl + 1):combined_nl, ]

    dos.F1_neutral <- genos.F1_neutral[, , 1] + genos.F1_neutral[, , 2]
    bed.F1_neutral <- as.bed.matrix(dos.F1_neutral)

    dos.F1_quanti <- genos.F1_quanti[, , 1] + genos.F1_quanti[, , 2]
    bed.F1_quanti <- as.bed.matrix(dos.F1_quanti)

    M<-matching(dos)
    k<-beta.dosage(M,MATCHING=TRUE)
    fst.founders<-fs.dosage(M,pop=rep(1:np,each=(ns+nd)),matching=TRUE)

    #matching for F1
    M.F1<-matching(dos.F1_neutral)

    #kinship estimated from markers for F1
    k.F1.g<-beta.dosage(M.F1,MATCHING =TRUE)
    #mean kinship of founders in each pop
    mean.pop.kinship<-unlist(lapply(0:(np-1),function(x) mean(mat2vec(k[x*ni+1:ni,x*ni+1:ni]))))
    nitF1<-nrow(k.F1.g)
    #Numbers of F1 per pop
    nipF1<-ns*nd*no
    k.F1.g.n<-matrix(0,nrow=nitF1,ncol=nitF1)
    # substract mean kinship per pop to inds in the pop, 
    # first founders and then F1
    for (i in 0:(np-1)){
        #founders
        k.F1.g.n[i*ni+1:ni,i*ni+1:ni]<-
        (k.F1.g[i*ni+1:ni,i*ni+1:ni]-mean.pop.kinship[i+1])/
        (1-mean.pop.kinship[i+1])
        #F1
        k.F1.g.n[nf,nf][i*nipF1+1:nipF1,i*nipF1+1:nipF1]<-
        (k.F1.g[nf,nf][i*nipF1+1:nipF1,i*nipF1+1:nipF1]-mean.pop.kinship[i+1])/
        (1-mean.pop.kinship[i+1])  
    }
    minFij<-min(mat2vec(fst.founders$FsM))
    Theta.P<-(fst.founders$FsM-minFij)/(1-minFij)

    The.M<-matrix(0,nrow=np*nipF1,ncol=np*nipF1)
    for (i in 0:(np-1)){
        The.M[i*nipF1+1:nipF1,i*nipF1+1:nipF1]<-
        hierfstat::kinship2grm(k.F1.g.n)[nf,nf][i*nipF1+1:nipF1,i*nipF1+1:nipF1]*(1-Theta.P[i+1,i+1])
    
    }
#####
    tp <- dim(The.M)[1]


    data_rep <- list(bed.F1_neutral = bed.F1_neutral,
                    tp = tp,
                    The.M = The.M,
                    Theta.P = Theta.P,
                    bed.F1_quanti = bed.F1_quanti,
                    dos.F1_quanti = dos.F1_quanti)

}

make.traits<-function(bed=bed,n.causal=1000,h2=0.5,minp=0.01,sel=0.0,pop=NULL,popid=NULL){

  ef.size <- rep(0.1,n.causal)
  gval<-as.numeric(gaston::as.matrix(bed)%*%ef.size)
  va<-var(gval)
  phenog<-gval+stats::rnorm(nrow(bed),sd=(va*(1/h2-1))^.5)
  if(!is.null(pop)){
    if (is.null(popid)) popid<-bed@ped$famid
    popid<-factor(popid)
    npop<-nlevels(popid)
    pop.efsize<-stats::rnorm(npop,sd=sd(phenog)*pop^0.5)
    for (i in 1:npop){
      x<-levels(popid)[i]
      phenog[popid==x]<-phenog[popid==x]+pop.efsize[i]
    }
  }
  phenog<-phenog-mean(phenog)
  h2hat<-va/var(phenog) #h2
  return(list(h2=h2,h2hat=h2hat,
              trait=data.frame(gval=gval-mean(gval),pheno=phenog),
              causal=data.frame(chr=bed@snps$chr,id=bed@snps$id,efs=ef.size)))
}

get_vars <- function(bed.F1_quanti, Theta.P, The.M, bed.F1_neutral, dos.F1_quanti) {
    total_f1<-length(bed.F1_neutral@ped$id[nf])
   
    ####
    
    #h2 <- 1
    #bed <- bed.F1_quanti[nf, ]
    #sampled_neutral <- sample(dim(bed.F1_neutral[nf,])[2],100)
    #bed_nm <- bed.F1_neutral[nf,sampled_neutral]
    #pheno <- make.traits(bed_nm, n.causal = 100, h2 = h2, minp = 0.001)
    #pheno <- make.traits(bed.F1_quanti[nf, ], n.causal = 100, h2 = h2, minp = 0.001)
    #Y <- pheno$trait$pheno
    ####Dos sum 
    dosage_quanti <- dos.F1_quanti[1:total_f1,] 
    Y <- rowSums((dosage_quanti-1)*0.02)
    err <- rnorm(length(Y), mean = 0, sd = 1)
    Y <- Y + err
    mean_Y <- mean(Y)
    centered_Y <- Y - mean_Y
    var_Y <- var(centered_Y)
    Y <- centered_Y / sqrt(var_Y)
    row.names(Theta.P) <- colnames(Theta.P) <- paste("pop", 1:np, sep = "_")
    row.names(The.M) <- colnames(The.M) <- bed.F1_neutral@ped$id[nf]
    two.Theta.P <- 2 * Theta.P
    ind_per_pop = length(bed.F1_neutral@ped$id[nf])/np
    dat <- data.frame(pop = paste("pop", rep(1:np, each = ind_per_pop), sep = "_"), ind = bed.F1_neutral@ped$id[nf], Y = Y)
    brms_mf <- brm(Y ~ 1 + (1 | gr(pop, cov = two.Theta.P)) + (1 | gr(ind, cov = The.M)), 
                   data = dat, data2 = list(two.Theta.P = two.Theta.P, The.M = The.M), 
                   family = gaussian(), chains = 8, cores = 4, iter = 3000, warmup = 1000, thin = 2)

    tmp <- data.frame(matrix(unlist(lapply(VarCorr(brms_mf, summary = FALSE),
                                           function(x) x$sd^2)), ncol = 3))
    names(tmp) <- names(VarCorr(brms_mf, summary = FALSE))

    quant_med <- quantile(tmp$pop - tmp$ind, c(0.5, 0.025, 0.975))
    mean_diff <- mean(tmp$pop - tmp$ind)
    
    hyp <- "sd_pop__Intercept^2 - sd_ind__Intercept^2 = 0"
    the_hyp <- hypothesis(brms_mf, hyp, class = NULL)

    post_samples <- posterior_samples(brms_mf, pars = c("sd_pop__Intercept", "sd_ind__Intercept"))
    post_samples$log_ratio <- log(post_samples$sd_pop__Intercept^2 / post_samples$sd_ind__Intercept^2)
    mean_log_ratio <- mean(post_samples$log_ratio)
    quant_log_ratio <- quantile(post_samples$log_ratio, probs = c(0.025, 0.975))

    
    
    frequency_opp <- 2 * sum(sign(post_samples$log_ratio) != sign(median(post_samples$log_ratio))) / length(post_samples$log_ratio)
    p_value <- 2 * min(frequency_opp, 1 - frequency_opp)


    results <- list(p_value = p_value,
                    BRMS_test = the_hyp$hypothesis[2:5],
                    BRMS_mean_diff = mean_diff,
                    BRMS.test.median = quant_med["50%"],
                    BRMS.test.median_lower = quant_med["2.5%"],
                    BRMS.test.median_upper = quant_med["97.5%"],
                    BRMS_mean_log_ratio = mean_log_ratio,
                    BRMS_log_ratio_ci_lower = quant_log_ratio[1],
                    BRMS_log_ratio_ci_upper = quant_log_ratio[2])
    return(results)
}

#Initializing lists to store our results
mean_list <- list()
ci_list <- list()
median_list <- list()
ci_median_list <- list()
log_list <- list()
ci_log_list <- list()
########################################


analyze_replicate <- function(replicate_number) {
    print(paste0("Processing replicate ", replicate_number))
    neutral_file <- paste0(Simulation_directory, Neutral_file_name, sprintf("%03d", replicate_number), ".dat")
    quanti_file <- paste0(Simulation_directory, Quanti_file_name, sprintf("%03d", replicate_number), ".dat")
    dat <- read.fstat(fname = neutral_file)
    dat_quanti <- read.fstat(fname = quanti_file)
    #dat <- subsampind(dat,sampsize=10)
    #dat_quanti <- subsampind(dat_quanti,sampsize=10)
    
    #Getting data for given replicate
    data_rep <- get_data(dat, dat_quanti)
    bed.F1_neutral <- data_rep$bed.F1_neutral
    tp <- data_rep$tp
    The.M <- data_rep$The.M
    Theta.P <- data_rep$Theta.P
    bed.F1_quanti <- data_rep$bed.F1_quanti
    dos.F1_quanti <- data_rep$dos.F1_quanti
    
    if (!all(eigen(The.M)$values > 0)) {
        The.M <- nearPD(The.M)$mat
        if (!all(eigen(The.M)$values > 0)){
            return(NULL)
        }
    }

    
    #Calculations for given replicate
    results <- get_vars(bed.F1_quanti, Theta.P, The.M, bed.F1_neutral, dos.F1_quanti) #BRMS 
    
    return(results)
}


#Analyze single replicate
results <- analyze_replicate(replicate_number)

#Append results to CSV
if (!is.null(results)) {
    results_df <- data.frame(replicate = replicate_number,
                             p_value = results$p_value,
                             mean_diff = results$BRMS_mean_diff,
                             CI_lower = results$BRMS_test$CI.Lower,
                             CI_upper = results$BRMS_test$CI.Upper,
                             median_diff = results$BRMS.test.median,
                             CI_median_lower = results$BRMS.test.median_lower,
                             CI_median_upper = results$BRMS.test.median_upper,
                             log_ratio = results$BRMS_mean_log_ratio,
                             CI_log_lower = results$BRMS_log_ratio_ci_lower,
                             CI_log_upper = results$BRMS_log_ratio_ci_upper)
    
    write.table(results_df, file = paste0(file_csv_name, ".csv"), append = TRUE, sep = ",", col.names = !file.exists(paste0(file_csv_name,".csv")), row.names = FALSE)
}
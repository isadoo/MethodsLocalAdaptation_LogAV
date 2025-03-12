.libPaths("/work/FAC/FBM/DEE/jgoudet/default/R/lib")
library(tidyverse)
library(hierfstat)
library(gaston)
library(JGTeach)
library(brms)
library(Matrix)
source("/users/ijeronim/mywork/Chapter1/driftsel.r")
# Set parameters
Population_structure <- "IM"
Selective_or_Neutral <- "Neutral"
replicate_number <- as.integer(commandArgs(trailingOnly = TRUE))
#replicate_number <- 3
generations <- 500
optima_type <- ""

# File paths
#Simulation_directory <- paste0("/work/FAC/FBM/DEE/jgoudet/pop_fst/isa/Chapter1/QN_SS_",Selective_or_Neutral,"_ALLFolders/QN_", Population_structure, "_", Selective_or_Neutral, optima_type, "_rep")
Simulation_directory <- paste0("/work/FAC/FBM/DEE/jgoudet/pop_fst/isa/Chapter1/QN_", Population_structure, "_", Selective_or_Neutral, optima_type, "/")
#Neutral_file_name <- paste0("/neutral_data_g", generations)
#Quanti_file_name <- paste0("/quanti_trait_g", generations)
Neutral_file_name <- paste0("neutral_data_g", generations, "_r")
Quanti_file_name <- paste0("quanti_trait_g", generations, "_r")

neutral_file <- paste0(Simulation_directory, Neutral_file_name, sprintf("%03d", replicate_number), ".dat")
quanti_file <- paste0(Simulation_directory, Quanti_file_name, sprintf("%03d", replicate_number), ".dat")
#neutral_file <- paste0(Simulation_directory,replicate_number, Neutral_file_name,".dat")
#quanti_file <- paste0(Simulation_directory, replicate_number,Quanti_file_name, ".dat")

# Naming output CSV
file_csv_name <- paste0("ComparisonAllMethods_Feb_", Population_structure, "_", Selective_or_Neutral, ".csv")

# Read data
sim <- read.fstat(fname = neutral_file)
sim_quanti <- read.fstat(fname = quanti_file)

# Process dosage data
dos<-biall2dos(sim[,-1])
dos_quanti <- biall2dos(sim_quanti[,-1], diploid = TRUE)

# FST for neutral markers
fst_neutral <- fs.dosage(dos, pop = sim[,1])$FsM

# Generate F1 generation and phenotype
# This part is reused for both methods
nqtl <- 100
ns <- 5
nd <- 5
np <- 8
no <- 2
n_ind <- no * ns * nd * np
nft <- np * (ns + nd)


# Functions used
create_pedigree <- function(nft, ns, nd, np, no) { 
  sire_F1 <- rep(1:ns, each = nd * no)
  dam_F1 <- rep(ns + 1:nd, each = no, ns)
  sire_F1 <- rep(0:(np - 1) * (nd + ns), each = (nd * ns * no)) + sire_F1
  dam_F1 <- rep(0:(np - 1) * (nd + ns), each = (nd * ns * no)) + dam_F1
  sire <- c(rep(NA, nft), sire_F1)
  dam <- c(rep(NA, nft), dam_F1)
  nt <- length(sire)
  nf <- -c(1:nft)
  ni <- ns + nd
  pop_P <- rep(1:np, each = (ns + nd))
  pop_F1 <- rep(1:np, each = (ns * nd * no))
  ped_forBreeding <- data.frame(ind = 1:nt, sire = sire, dam = dam)
  ped <- data.frame(id = (max(dam_F1) +1):((ns*nd*no*np)+max(dam_F1)), sire = sire_F1 , dam = dam_F1, sire.pop = pop_F1, dam.pop = pop_F1)
  pedigrees <- list(ped = ped,ped_forBreeding = ped_forBreeding)
}



generate_f1_combined <- function(sim, sim_quanti, ped_forBreeding) {
  dos<-biall2dos(sim[,-1])
  nl<-ncol(dos)

      #dos_quanti <- fstat2dos(dat_quanti[, -1])
  dos_quanti <- biall2dos(sim_quanti[, -1])
  nl_quanti <- ncol(dos_quanti)

  combined_nl <- nl + nl_quanti
  combined_founders <- array(0, dim = c(nft, combined_nl, 2))
  combined_dat <- cbind(sim[, -1], sim_quanti[, -1])

  tmp1<-as.matrix(combined_dat%/%10-1)
  tmp2<-as.matrix(combined_dat%%10-1)
  combined_founders[,,1]<-tmp1
  combined_founders[,,2]<-tmp2

  nfounders_combined <- rbind(combined_founders[, , 1], combined_founders[, , 2])

  genos.F1_combined <- drop.along.ped(ped_forBreeding, founders.genotypes = nfounders_combined, nloc = combined_nl, maplength = 1000)

  #Extract the neutral and quantitative genotypes from the combined result
  genos.F1_neutral <- genos.F1_combined[, 1:nl, ]
  genos.F1_quanti <- genos.F1_combined[, (nl + 1):combined_nl, ]

  dos.F1_neutral <- genos.F1_neutral[, , 1] + genos.F1_neutral[, , 2]
  bed.F1_neutral <- as.bed.matrix(dos.F1_neutral)

  dos.F1_quanti <- genos.F1_quanti[, , 1] + genos.F1_quanti[, , 2]
  bed.F1_quanti <- as.bed.matrix(dos.F1_quanti)
  
  list(
    neutral = dos.F1_neutral,
    quanti = dos.F1_quanti,
    neutral_bed = bed.F1_neutral, 
    quanti_bed = bed.F1_quanti
  )
}

create_data_frame <- function(n_ind, np, ns, nd, no, phenotype) {
  nft <- np * (ns + nd)
  individuals <- nft:(n_ind+nft-1)
  populations <- rep(1:np, each = ns * nd * no)
  sires <- rep(rep(1:ns, each = nd * no), np)
  dams <- rep(rep(1:nd, each = no), ns * np)
  dataf <- data.frame(Individual = individuals, Population = factor(populations), Sire = factor(sires), Dam = factor(dams), Y = phenotype)
}

run_qstfst <- function(dataf, fst_neutral, ns, nd, no, np) {
  nrep = 1000
  popn<-factor(rep(1:np,each=(ns*nd*no)))
  siren<-rep(1:ns,each=nd*no)
  damn<-rep(1:nd,each=no,ns)

  #y<-rnorm(ns*nd*no*np)
  #dat<-data.frame(Y=Y[nf],pop=popn,sire=factor(rep(siren,np)),dam=factor(rep(damn,np)))
  my.anova<-anova(aov(Y~Population/(Sire+Dam),data=dataf))  #assuming no dominance
  MSs<-my.anova[,3]
  DFs<-my.anova[,1]
  MSsire<-MSs[2]
  MSdam<-MSs[3]
  MSpop<-MSs[1]
  MSwithin<-MSs[4]
  DFsire<-DFs[2]
  DFdam<-DFs[3]
  DFpop<-DFs[1]
  DFwithin<-DFs[4]

  #to modify if needs be
  Fst<-fst.founders$Fs[2,np+1]

  get.rand.fst<-function(dos,pop){
    nl<-ncol(dos)
    np<-length(table(pop))
    x<-sample(nl,size=nl,replace=TRUE)
    hierfstat::fst.dosage(dos[,x],pop=pop)[np+1]
    }

  fst.hat<-replicate(nrep,get.rand.fst(dos,rep(1:np,each=ni)))

  #to get MSPopNeutral
  x.hat<-fst.hat/(1-fst.hat)

  MSsire.hat<-MSsire/DFsire*rchisq(nrep,DFsire)
  MSdam.hat<-MSdam/DFdam*rchisq(nrep,DFdam)
  MSwithin.hat<-MSwithin/DFwithin*rchisq(nrep,DFwithin)

  #estimated variance components for dams and sires
  sigA.dam<-(MSdam-MSwithin)/(ns*no)
  sigA.sire<-(MSsire-MSwithin)/(nd*no)

  sigA.dam.hat<-(MSdam.hat-MSwithin.hat)/(ns*no)
  sigA.sire.hat<-(MSsire.hat-MSwithin.hat)/(nd*no)

  #needed?
  sigA.hat<-(sigA.dam.hat+sigA.sire.hat)/2
  VA.hat<-4*sigA.hat

  #only place to get variablity for fst I think
  MSpopNeutral<-MSwithin+sigA.dam*ns*no+sigA.sire*nd*no+4*sigA.dam*(ns*nd*no*Fst/(1-Fst))+4*sigA.sire*(nd*ns*no*Fst/(1-Fst))
  MSpopneutral.hat<-MSpopNeutral/DFpop*rchisq(nrep,df=DFpop)  

  #estimated component of variance for pop
  sigA.pop<-(MSpop-sigA.dam*no*ns-sigA.sire*no*nd-MSwithin)/(no*ns*nd)

  #estimated qst, assuming no maternal nor dominance effects
  QST<-sigA.pop/(sigA.pop+4*sigA.sire+4*sigA.dam)


  sigA.pop.neutral.hat<-(MSpopneutral.hat-sigA.dam.hat*no*ns-sigA.sire.hat*no*nd-MSwithin.hat)/(no*ns*nd)

  #the null distribution for QST
  #could be easily expanded at no computing costs
  #by randomizing independently all the sources of variations

  QST.neutral.hat<-sigA.pop.neutral.hat/(sigA.pop.neutral.hat+4*sigA.dam.hat+4*sigA.sire.hat)

  tmp.qstwg<-data.frame(fst.star=fst.hat,qstneut.star=QST.neutral.hat,sApopneut.star=sigA.pop.neutral.hat,
  sASire.star=sigA.sire.hat,sAdam.star=sigA.dam.hat,sAwithin.star=MSwithin.hat)


  five.num.qstgw<-quantile(QST.neutral.hat-QST- fst.hat + Fst,c(0.025,0.25,0.5,0.75,0.975))
  pval.qstgw.neg<-(sum(QST.neutral.hat-fst.hat <= QST- Fst )+1)/(nrep+1)
  pval.qstgw.pos<-(sum(QST.neutral.hat - fst.hat >= QST - Fst)+1)/(nrep+1)

  p_value <- sum(abs(QST.neutral.hat-fst.hat) >= abs(QST- Fst)) / length(QST.neutral.hat-fst.hat)
  #combined_p_values <- c(pval.qstgw.neg, pval.qstgw.pos)
  return(data.frame(method = "QST-FST", p_value = p_value))
}

run_lava <- function(dataf,Theta.P,The.M) {
#####
  #tp <- dim(The.M)[1]

  row.names(Theta.P) <- colnames(Theta.P) <- paste("pop", 1:np, sep = "_")
  row.names(The.M) <- colnames(The.M) <- dataf$Individual
  two.Theta.P <- 2 * Theta.P
  ind_per_pop = length(dataf$Individual)/np
  dat <- data.frame(pop = paste("pop", rep(1:np, each = ind_per_pop), sep = "_"), ind = dataf$Individual, Y = dataf$Y)
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

    
    
  p_value <- 2 * sum(sign(post_samples$log_ratio) != sign(median(post_samples$log_ratio))) / length(post_samples$log_ratio)
    
  
  return(data.frame(method = "LAVA", p_value = p_value, log_ratio = mean_log_ratio))

}

run_driftsel <- function(dataf,dos,ped,covars){
  boot_coan <- function(list_blocks,bed,np,indperpop){
    sampled_blocks <- sample(list_blocks, 100, replace = TRUE) # default 100 = number of blocks
    samp_list<-unlist(sampled_blocks)
    fst_mat <- fs.dosage(bed[,samp_list],pop=rep(1:np,each=indperpop))$FsM
    off<- fst_mat[row(fst_mat)!=col(fst_mat)]
    coan <- (fst_mat - min(off))/(1 - min(off))
  }
  
  bed <- as.bed.matrix(dos)
  n_snps <- length(bed@p)

  block <- 1:100
  num_blocks <- floor(n_snps/100)
  list_blocks <- lapply(1:num_blocks, function(i) {
    current_block <- block + 100 * (i - 1)
  })

  tensor.Theta.P <- replicate(25, boot_coan(list_blocks,bed, np, indperpop = ns+nd)) #25 is the number of bootstraps
  traits <- data.frame(id = 1:(ns*nd*no*np), phenotype = dataf$Y)
  
  samp <- MH(tensor.Theta.P, ped, covars, traits,10000, 5000, 20, alt=T)
  s <- neut.test(samp$pop.ef, samp$G, samp$theta,, silent=F)
  driftsel_results <- list(S_Value = mean(s))
}

Get_TheM_ThetaP <- function(dos, dos.F1_neutral){
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


  #if (!all(eigen(The.M)$values > 0)) {
  #      The.M <- nearPD(The.M)$mat
  #      if (!all(eigen(The.M)$values > 0)){
  #          return(NULL)
  #      }
  #  }

  #if (!all(eigen(Theta.P)$values > 0)) {
  #      Theta.P <- nearPD(Theta.P)$mat
  #      if (!all(eigen(Theta.P)$values > 0)){
  #          return(NULL)
  #      }
  #  }

  coancestries <- list(The.M = The.M, Theta.P = Theta.P)
}

#Pedigree
ped_list <- create_pedigree(nft, ns, nd, np, no)
ped <- ped_list$ped
ped_forBreeding <- ped_list$ped_forBreeding
nt<-length(ped$sire)
nf<- -c(1:nft)
ni<-ns+nd

#Generate F1 genotypes
combined_genos <- generate_f1_combined(sim, sim_quanti, ped_forBreeding)

#Extract neutral and quantitative genotypes
dos.F1_neutral <- combined_genos$neutral
dos.F1_quanti <- combined_genos$quanti
bed.F1_neutral <- combined_genos$neutral_bed
bed.F1_quanti <- combined_genos$quanti_bed

dos.Just.F1_quanti <- dos.F1_quanti[(nft+1):dim(dos.F1_quanti)[1],]
phenotype = rowSums((dos.Just.F1_quanti-1)*0.2)+rnorm(nrow(dos.Just.F1_quanti))

#Create dataframe for analysis
dataf <- create_data_frame(n_ind, np, ns, nd, no, phenotype)

#Run QST-FST method
pop_P <- rep(1:np, each = (ns + nd))
fst.founders <- fs.dosage(dos,pop_P)
qstfst_result <- run_qstfst(dataf, fst.founders, ns, nd, no, np)

#Get M and Theta P coancestry matrices
coancestries <- Get_TheM_ThetaP(dos, dos.F1_neutral)
The.M <- coancestries$The.M
Theta.P <- coancestries$Theta.P
#Run LAVA method
lava_result <- run_lava(dataf,Theta.P, The.M)

#Run Driftsel method
covars <- data.frame(id = 1:(ns*nd*no*np), dummy = rep(1,(ns*nd*no)))
driftsel_result <- run_driftsel(dataf,dos,ped,covars)

#Combine results
results_df <- data.frame(
  replicate_number = replicate_number,
  p_value_QSTFST = qstfst_result$p_value,
  p_value_LAVA = lava_result$p_value,
  log_ratio_LAVA = lava_result$log_ratio,
  S_value_Driftsel = driftsel_result$S_Value
)

write.table(
  results_df, file = file_csv_name, append = TRUE, 
  sep = ",", col.names = !file.exists(file_csv_name), row.names = FALSE
)
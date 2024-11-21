
.libPaths("/work/FAC/FBM/DEE/jgoudet/default/R/lib")
library(tidyverse)
library(hierfstat)
library(boot)
library(JGTeach)
library(gaston)
library(parallel)
library(MCMCpack)

source("driftsel.r")

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
file_csv_name <- paste0("MethodsLocalAdaptation_DriftselTT_", Selective_or_Neutral, "_", Population_structure)

#Default breeding parameters
ns <- 5  # Number of sires
nd <- 5  # Number of dams
no <- 2  # Number of offspring per pair
n_replicates <- 100  # Number of replicates
n_loci_qtl <- 100  # Number of QTLs
maplength_quanti <- 50  # Map length for drop along ped
##################################################################################################

#Simulation data--------------------------------------------------

neutral_file <- paste0(Simulation_directory, Neutral_file_name, sprintf("%03d", replicate_number), ".dat")
quanti_file <- paste0(Simulation_directory, Quanti_file_name, sprintf("%03d", replicate_number), ".dat")
sim <- read.fstat(fname = neutral_file)
sim_quanti <- read.fstat(fname = quanti_file)
dos <- biall2dos(sim[,-1],diploid=TRUE)
dos_quanti <- biall2dos(sim_quanti[,-1],diploid=TRUE)
FsM<-fs.dosage(dos,pop=sim[,1])$FsM #
n_ind <- no * ns * nd * np  #total number of individuals F1
indperpop <- 1000 #Original population size (before sampling)
ThetaTensor_rep <- 25
#-------------------------------------------------------


#Pedigree-------
nft <- np * (ns + nd)
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

#-----------------
#FST for later use
fst.founders <- fs.dosage(dos,pop_P)

#Building F1 generation ----------------------------------------------------------

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
#------------------------------------------------------------------------

#Building phenotype
dos.Just.F1_quanti <- dos.F1_quanti[(nft+1):dim(dos.F1_quanti)[1],]
phenotype = rowSums((dos.Just.F1_quanti-1)*0.2)+rnorm(nrow(dos.Just.F1_quanti))

#Building phenotype
#dos.Just.F1_quanti <- dos.F1_quanti[(nft+1):dim(dos.F1_quanti)[1],]
#phenotypes = rowSums((dos.Just.F1_quanti-1)*0.2)+rnorm(nrow(dos.Just.F1_quanti))
#traits <- data.frame(id = 1:(ns*nd*no*np), trait.1 = phenotypes)
#dos.Just.F1_neutral <- dos.F1_neutral[(nft+1):dim(dos.F1_neutral)[1],]
covars <- data.frame(id = 1:(ns*nd*no*np), dummy = rep(1,(ns*nd*no)))

#------------Building Theta
MatchingDos <- matching(dos)
fst.founders<-fs.dosage(MatchingDos,pop=rep(1:np,each=(ns+nd)),matching=TRUE)
minFij<-min(mat2vec(fst.founders$FsM))
Theta.P<-(fst.founders$FsM-minFij)/(1-minFij)
#From matrix to tensor
#Same as afm$Theta = tensor.Theta.P from driftsel manual

#Function to bootstrap Theta.P
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


tensor.Theta.P <- replicate(ThetaTensor_rep, boot_coan(list_blocks,bed, np, indperpop = ns+nd))

s_means <- c()
n_traits = 100
nc = 100

total_F1<-dim(dos.Just.F1_quanti)[1]
traits <- data.frame(id = 1:(ns*nd*no*np), phenotype = phenotype)

#DRIFTSEL FUNCTIONS PART ################################################
samp <- MH(tensor.Theta.P, ped, covars, traits,10000, 5000, 20, alt=T)
s <- neut.test(samp$pop.ef, samp$G, samp$theta,, silent=F)
##########################################################################
results_df <- data.frame(ReplicateNumber = replicate_number, S_Value = s)
file_csv_name <- paste0("DRIFTSEL_QuantiData_",Population_structure)


#dataframe to CSV
write.table(results_df, file = paste0(file_csv_name, ".csv"), append = TRUE, sep = ",", col.names = !file.exists(paste0(file_csv_name,".csv")), row.names = FALSE)

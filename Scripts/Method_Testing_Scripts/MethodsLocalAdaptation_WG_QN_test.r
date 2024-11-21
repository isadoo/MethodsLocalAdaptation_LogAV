.libPaths("/work/FAC/FBM/DEE/jgoudet/default/R/lib")
library(tidyverse)
library(hierfstat)
library(boot)
library(JGTeach)
library(gaston)


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
file_csv_name <- paste0("MethodsLocalAdaptation_QSTFST_", Selective_or_Neutral, "_", Population_structure)

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
#-------------------------------------------------------


#Pedigree-------
nft <- np * (ns + nd)
sire <- rep(1:ns, each = nd * no)
dam <- rep(ns + 1:nd, each = no, ns)
sire <- rep(0:(np - 1) * (nd + ns), each = (nd * ns * no)) + sire
dam <- rep(0:(np - 1) * (nd + ns), each = (nd * ns * no)) + dam
sire <- c(rep(NA, nft), sire)
dam <- c(rep(NA, nft), dam)
nt <- length(sire)
nf <- -c(1:nft)
ni <- ns + nd
pop_P <- rep(1:np, each = (ns + nd))
pop_F1 <- rep(1:np, each = (ns * nd * no))
ped <- data.frame(ind = 1:nt, sire = sire, dam = dam)
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

genos.F1_combined <- drop.along.ped(ped, founders.genotypes = nfounders_combined, nloc = combined_nl, maplength = 1000)

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

#DATA
individuals <- 1:n_ind
populations <- rep(1:np, each = ns * nd * no)
sires <- rep(rep(1:ns, each = nd * no), np)
dams <- rep(rep(1:nd, each = no), ns * np)

dataf <- data.frame(Individual = individuals, Population = factor(populations), Sire = factor(sires), Dam = factor(dams), Y = phenotype)

# --------------------------------------------------------------------------------------


#JEROME'S PORTION OF THE SCRIPT: 
nboot <- 1000
nrep<-nboot

##setting up the anova for NCII
##calculates Guillaume-Gilbert-Whitlock QST-FST test
##adapted to NCII design assuming no sireXdam effects


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

results_df <- data.frame(replicate_number = replicate_number,
                         estimated_QST = QST,
                         p_value_pos = pval.qstgw.pos,
                         p_value_neg = pval.qstgw.neg,
                         p_value = p_value)

write.table(results_df, file = file_csv_name, append = TRUE, sep = ",", col.names = !file.exists(file_csv_name), row.names = FALSE)

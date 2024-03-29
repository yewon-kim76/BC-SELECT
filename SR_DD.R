### Breast Cancer TCGA + METABRIC + updated cox regression filters
### Molecular Taxonomy of Breast Cancer International Consortium (METABRIC)
### TCGA Immune Checkpoint Inhibitor (TCGA.ICI)
### BC TCGA + METABRIC training
### Remove SCNA/CTLA4p, update PD1/PDL1

# Training---------------------------------------------------------------------------------------------------
rm(list=ls())
library(readxl); library(rmarkdown);library(ggplot2); library(tidyverse); library(dplyr); library(pROC); library(MKmisc)
library(hrbrthemes); library(viridis); library(forcats); library(data.table); library(Rcpp); 
library(survival); library(parallel); library(pracma); library(ROCR);
set.seed(12345)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")
setwd("/gpfs/gsfs12/users/kimy35/SELECT/select_targeted/cox_sort_12_05/R_immuno_data")
source("./source2.R")

# Read data-----------------------------------------------------------------------------------------
metabric.list <- readRDS("./metabric.list.rds") # METABRIC RNA-seq
metabric.surv <- readRDS("./metabric.surv.rds") # METABRIC survival phenotypic data
TCGA.list <- readRDS("./TCGA.list.rds") # TCGA RNA-seq
TCGA.surv <- readRDS("./TCGA.surv.rds") # TCGA survival phynotypic data
RNA_data <- readRDS("./RNA_data.rds") # Aggregated RNA-seq (TCGA + METABRIC)

nanostring=fread("/gpfs/gsfs12/users/kimy35/SELECT/select_immuno/data/nanostring.txt")$"Gene Name" ## 1154 genes

targets = c("PD1/PDL1") 
target.genes = c("CD274","PDCD1") # Both of them used in phylogenetic profiling
targetIDs = match(targets, RNA_data$genes) # 19002
partnerIDs = match(nanostring, RNA_data$genes) # 58 missing
partnerIDs = partnerIDs[!is.na(partnerIDs)] # Remove genes not in TCGA gene list
partner.genes = RNA_data$genes[partnerIDs] # Gene name corresponding to partnerIDs
sr = cbind(partnerIDs,targetIDs) 
sr.genes = cbind(rep(partner.genes,2), rep(target.genes, each=1096))## Phylogenetic

# (1-1) Tumor screen (hypergeometric): combined mRNA 
sr.curr = sr ##partnerID, targetID
sourceCpp("/gpfs/gsfs12/users/kimy35/SELECT/select_immuno/R/HyperGeometricTest.pair.cpp",rebuild=T)
pval.mRNA1 = hypergeometricTestPair(scnaq= RNA_data$mRNAq2, pairs=sr.curr) 
pval.mRNA.up = hypergeometricTestPair(scnaq= RNA_data$mRNAq2, pairs=sr.curr, lowerTail=0)
pval.mRNA = cbind(sr.curr[,1:2], pval.mRNA1, pval.mRNA.up)
colnames(pval.mRNA)[c(3:11)] <- c("pval.mRNA1_1","pval.mRNA1_2","pval.mRNA1_3","pval.mRNA1_4",
                                  "pval.mRNA1_5","pval.mRNA1_6","pval.mRNA1_7","pval.mRNA1_8",
                                  "pval.mRNA1_9")
colnames(pval.mRNA)[c(12:20)] <- c("pval.mRAN.up_1","pval.mRAN.up_2","pval.mRAN.up_3",
                                   "pval.mRAN.up_4",  "pval.mRAN.up_5","pval.mRAN.up_6","pval.mRAN.up_7",
                                   "pval.mRAN.up_8","pval.mRAN.up_9")

#plot(pval.mRNA[,3],pval.mRNA[,12]) ## 1-pval.mRNA[,3] = pval.mRNA[,12] perfect negative correlation
rm(pval.mRNA1, pval.mRNA.up)

# (1-2) Tumor screen (cox regression): for each dataset
coxph.run = function(formula.list, data1){
  fmla <- as.formula(paste(" Surv(time,status) ~ ", paste(formula.list, collapse= "+")))
  cox.out = coxph(fmla, data=data1); 
  aa  = summary(cox.out); out = aa$coefficients["cov", ] ## saved values = "cov"
  out}

coxph.robust = function(data1, f1.list) {coxph.run(c(f1.list), data1)} 

sr.clinical.screen = function(pair, prob, pheno, f1.list=NULL){
  g1 = prob$mRNA.norm[pair[1],] # partner ID 
  g2 = prob$mRNA.norm[pair[2],] # target ID (PD1/PDL1)
  f1 = prob$mRNAq2[pair[1],]  # partner ID
  f2 = prob$mRNAq2[pair[2],]  # target ID (PD1/PDL1)
  surv.strata = data.frame(pheno$time, pheno$status); colnames(surv.strata) <- c("time","status")
  
  var <- pheno %>% select(any_of(c("types","age","sex","race","stage")))
  f1.list = names(var) # Select any possible variables
  f1.list <- ifelse(f1.list %in% c("types","sex","race","stage"), paste("strata(",f1.list,")", sep="" ), f1.list) # Change to strata
  surv.strata = data.frame(surv.strata, var)
  
  dt1 = cbind(surv.strata, cbind(g1 , g2))
  cov = ifelse(f2 == 0 & f1 == 0, 1, 0 ) # SR_DD
  dt1$cov = qnorm.array(cov)
  cox.out2 = coxph.robust(f1.list = c("g1",  "g2",  f1.list, "cov"),  data1=dt1) 
  uu = c(cox.out2)
  return(uu)} 

# TCGA cox regression with 4 variables
cox.mRNA = mclapply(1:nrow(sr.curr), function(tt) { 
  rs=rep(NA,10)    
  sr.tmp=sr.curr[tt,]; prob = TCGA.list # check the input
  if(sum(prob$mRNAq2[sr.tmp[2],]==0 & prob$mRNAq2[sr.tmp[1],]==0, na.rm=T)>=1 )
    rs=sr.clinical.screen(sr.curr[tt,], prob = TCGA.list, pheno = TCGA.surv) ##
  return(rs)}, mc.cores=32)
cox.du.mRNA.curr = do.call(rbind, cox.mRNA)
cox.du.mRNA.curr = cox.du.mRNA.curr[,1:5]
colnames(cox.du.mRNA.curr) <- c("bb_coef","bb_exp(coef)","bb_se(coef)","bb_z","bb_pvalue")
cox.du.mRNA=cbind(sr.curr,cox.du.mRNA.curr)
rm(cox.du.mRNA.curr, cox.mRNA)

clinical.beta = cox.du.mRNA[,c(3)] # Coefficient for bb
class(clinical.beta) = "numeric"

clinical.p = cox.du.mRNA[,c(7)] # p-values for bb (un-adjusted p-value)
class(clinical.p) = "numeric"

TCGA_clinp=cbind(clinical.beta,clinical.p)
colnames(TCGA_clinp) <- c("bb_coef_tcga","bb_pvalue_tcga")
rm(clinical.beta, clinical.p, cox.du.mRNA)

# METABRIC cox regression with 2 variables
cox.mRNA = mclapply(1:nrow(sr.curr), function(tt) { 
  rs=rep(NA,10)    
  sr.tmp=sr.curr[tt,];prob = metabric.list # check the input
  if(sum(prob$mRNAq2[sr.tmp[2],]==0 & prob$mRNAq2[sr.tmp[1],]==0,na.rm=T)>=1 )
    rs=sr.clinical.screen(sr.curr[tt,], prob = metabric.list, pheno = metabric.surv) ##
  return(rs)}, mc.cores=32)
cox.du.mRNA.curr = do.call(rbind, cox.mRNA)
cox.du.mRNA.curr = cox.du.mRNA.curr[,1:5]
colnames(cox.du.mRNA.curr) <- c("bb_coef","bb_exp(coef)","bb_se(coef)","bb_z","bb_pvalue")
cox.du.mRNA=cbind(sr.curr,cox.du.mRNA.curr)
rm(cox.du.mRNA.curr, cox.mRNA)

clinical.beta = cox.du.mRNA[,c(3)] # Coefficient for bb
class(clinical.beta) = "numeric"

clinical.p = cox.du.mRNA[,c(7)] # p-values for bb (un-adjusted p-value)
class(clinical.p) = "numeric"

metabric_clinp=cbind(clinical.beta,clinical.p)
colnames(metabric_clinp) <- c("bb_coef_meta","bb_pvalue_meta")
rm(clinical.beta, clinical.p, cox.du.mRNA)

# (2) Phylogenetic screen
load("./phylogenetic.profile.RData")
load("./feature.weight.RData")

sr.genes = data.table(rescuer = sr.genes[,1], vulnerable = sr.genes[,2])
phylogenetic.distance.score = phylo.profile(sr.genes)
phylogenetic.distance.score = cbind(phylogenetic.distance.score[1:1096], 
                                    phylogenetic.distance.score[1097:2192])

# Save all trained information
cox_all <- cbind(TCGA_clinp, metabric_clinp)
sr.tot=cbind(sr[,c(2,1)],pval.mRNA[,c(4,5,12)], cox_all, phylogenetic.distance.score) 
colnames(sr.tot)[10] <- "phylogenetic.distance.score_CD274"
colnames(sr.tot)[11] <- "phylogenetic.distance.score_PDCD1"

#saveRDS(sr.tot, file = "SR_DD.rds")

# Evaluation---------------------------------------------------------------------------------------------------
# (3) Find SR partner gene pairs, make a decision/classification (prediction)
rm(list=ls())
library(readxl); library(rmarkdown);library(ggplot2); library(tidyverse); library(dplyr); library(pROC); library(MKmisc)
library(viridis); library(forcats); library(data.table); library(Rcpp); 
library(survival); library(parallel); library(pracma); library(ROCR);
set.seed(12345)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")
setwd("/gpfs/gsfs12/users/kimy35/SELECT/select_targeted/cox_sort_12_05/R_immuno_data")
source("./source2.R")

sr.tot <- readRDS("./SR_DD.rds")
genes <- readRDS(file = "./TCGA_genes.rds")
#bc_gene_list <- read.delim("./bc_gene_list.txt", stringsAsFactors=TRUE)
#bc_gene_list <- as.character(bc_gene_list$Gene.Symbol)
nanostring=fread("/gpfs/gsfs12/users/kimy35/SELECT/select_immuno/data/nanostring.txt")$"Gene Name" ## 1154 genes

YALE <- readRDS("./YALE.rds") # Durvalumab, NCT02489448
ISPY_pemb <- readRDS("./ISPY_pemb.rds") # Pembrolizumab, ISPY-2
GSE173839_ISPY2 <- readRDS("./GSE173839_ISPY2.rds") # Durvalumab, ISPY-2

##
sr.final=NULL; tgs=c("PD1/PDL1") # Current target: anti-PD1/PDL1
ix = which(as.numeric(sr.tot[,1]) %in% match(tgs, genes))
sr.tot1 = sr.find.dd(sr.tot[ix,]) # 10 genes pass 3 filters
rnk=sr.ranking.dd(sr.tot1) # Sort the genes based on min phylogetic scores (PDCD1, CD274)
thr = 10; sr.final=rbind(sr.final, sr.tot1[rnk[1:thr],]) 
genes[unname(sr.final[,2])]
# "CEP55" "PPARGC1B" "TMEM173" "LILRB3" "C2" "LAG3" "TTK" "APH1B" "CCND3" "REN" 

####### Original SELECT SR partner genes
SELECT <- c("CXCL16", "ICAM4","LTBR","CD8A","IFITM2",
            "CD27","CD4","IL15RA","TNFRSF13B","TNFRSF13C")

#genes[which(genes %in% SELECT)]
sr.final = sr.tot1 <- sr.tot[which(sr.tot[,2] %in% which(genes %in% SELECT)), ]
genes[sr.tot1[,2]]
####### BC-SELECT better than original SELECT SR partner genes

# NCT02489448 (YALE)
auc=lgp1=eff1=rep(NA,1); sco=NULL;sco$pval=NA
ires1 = which(YALE$res_yale$binary == "1") # index of response
iirs1 = which(YALE$res_yale$binary == "0") # index of non-response
sco = eval.auc.dd(sr.final,YALE,ires1,iirs1)
s_yale_AUC1 = getAUC(sco$score, sco$flag)[1] ; s_yale_sr_genes_10 <- sco$sr.x[,2]; print(round(s_yale_AUC1,3))
s_yale_AUC2 = getAUC(sco$score2, sco$flag)[1]; print(round(s_yale_AUC2,3))
s_yale_AUC3 = getAUC(sco$score3, sco$flag)[1]; print(round(s_yale_AUC3,3))

# Area under the PR curve
round(getAUC(sco$score, sco$flag)[2],3) # GI only
round(getAUC(sco$score2, sco$flag)[2],3) # BC-SELECT
round(getAUC(sco$score3, sco$flag)[2],3) # Target

# GSE194040 Pembrolizumab, ISPY-2 
auc=lgp1=eff1=rep(NA,1); sco=NULL;sco$pval=NA
ires1 = which(ISPY_pemb$res_ispy_pemb$pCR == "1") # index of response
iirs1 = which(ISPY_pemb$res_ispy_pemb$pCR == "0") # index of non-response
sco = eval.auc.dd(sr.final,ISPY_pemb,ires1,iirs1)
s_ispy_pemb_AUC1 = getAUC(sco$score, sco$flag)[1] ; s_ISPY_pemb_sr_genes_10 <- sco$sr.x[,2]; print(round(s_ispy_pemb_AUC1,3))
s_ispy_pemb_AUC2 = getAUC(sco$score2, sco$flag)[1]; print(round(s_ispy_pemb_AUC2,3))
s_ispy_pemb_AUC3 = getAUC(sco$score3, sco$flag)[1]; print(round(s_ispy_pemb_AUC3,3))

# Area under the PR curve
round(getAUC(sco$score, sco$flag)[2],3) # GI only
round(getAUC(sco$score2, sco$flag)[2],3) # BC-SELECT
round(getAUC(sco$score3, sco$flag)[2],3) # Target


# GSE173839, Durvalumab, ISPY-2
auc=lgp1=eff1=rep(NA,1); sco=NULL;sco$pval=NA
ires1 = which(GSE173839_ISPY2$res_GSE173839_patient$binary == "1") # index of response
iirs1 = which(GSE173839_ISPY2$res_GSE173839_patient$binary == "0") # index of non-response
sco = eval.auc.dd(sr.final,GSE173839_ISPY2,ires1,iirs1)
s_ispy_durv_AUC1 = getAUC(sco$score, sco$flag)[1]; s_ISPY_Durv_sr_genes_10 <- sco$sr.x[,2]; print(round(s_ispy_durv_AUC1,3))
s_ispy_durv_AUC2 = getAUC(sco$score2, sco$flag)[1]; print(round(s_ispy_durv_AUC2,3))
s_ispy_durv_AUC3 = getAUC(sco$score3, sco$flag)[1]; print(round(s_ispy_durv_AUC3,3))

# Area under the PR curve
round(getAUC(sco$score, sco$flag)[2],3) # GI only
round(getAUC(sco$score2, sco$flag)[2],3) # BC-SELECT
round(getAUC(sco$score3, sco$flag)[2],3) # Target


# AUC Null distribution for checking significance
# NCT02489448 (YALE)
yale_bcAUC1 <- NULL; yale_bcAUC2 <- NULL # Saved AUC
bc_yale_genes <- YALE$genes[-which(YALE$genes %in% s_yale_sr_genes_10)] # Remove SR partners
bc_yale_genes <- bc_yale_genes[! bc_yale_genes %in% c('PDCD1', 'CD274')] # Remove target gene exp

for(i in 1:5000){ 
  pair <- sample(bc_yale_genes, size = 10, replace = FALSE, prob = NULL) 
  null_sr_score <- null_sr.score(YALE, pair)
  pval <- null_sr_score[c(which(YALE$res_yale$binary == "1"),
                          which(YALE$res_yale$binary == "0"))] # Reorder
  flag <- c(rep(1,25),rep(0,25)) ; target <- target_value(YALE)
  target <- target[c(which(YALE$res_yale$binary == "1"),
                     which(YALE$res_yale$binary == "0"))] # Reorder
  yale_bcAUC1[i] <- getAUC(pval, flag)[1] # GI only
  yale_bcAUC2[i] <- getAUC( (pval* rank.norm(target)), flag)[1] } # BC-SELECT

hist(yale_bcAUC1); hist(yale_bcAUC2)
round(pnorm(s_yale_AUC1, mean = mean(yale_bcAUC1), sd = sd(yale_bcAUC1), lower.tail=FALSE),4) # GI only
round(pnorm(s_yale_AUC2, mean = mean(yale_bcAUC2), sd = sd(yale_bcAUC2), lower.tail=FALSE),4) # BC-SELECT

# Empirical p-value: target gene expression case
control_bc <-NULL;
control_gene_bc <- YALE$genes[-which(YALE$genes %in% s_yale_sr_genes_10)] # Remove SR partners
control_gene_bc <- control_gene_bc[! control_gene_bc %in% c('PDCD1', 'CD274')] # Remove target gene exp

for(i in 1:5000){
  nn_gene <- sample(control_gene_bc, size = 2, replace = FALSE, prob = NULL) 
  pval_bc <- YALE$mRNA[nn_gene, c(which(YALE$res_yale$binary == "1"), which(YALE$res_yale$binary == "0"))] # Reorder
  pval_bc2 <- as.numeric(pval_bc[1,])*as.numeric(pval_bc[2,])
  flag <- c(rep(1,25),rep(0,25)); control_bc[i] <- getAUC(rank.norm(pval_bc2), flag)[1]
  rm(nn_gene)}

hist(control_bc) 
round(pnorm(s_yale_AUC3, mean = mean(control_bc), sd = sd(control_bc), lower.tail=FALSE),4) # Target


# GSE194040 Pembrolizumab, ISPY-2 
ISPY_bcAUC1 <- NULL; ISPY_bcAUC2 <- NULL # Saved AUC
bc_ISPY_genes <- ISPY_pemb$genes[-which(ISPY_pemb$genes %in% s_ISPY_pemb_sr_genes_10)] # Remove SR partners
bc_ISPY_genes <- bc_ISPY_genes[! bc_ISPY_genes %in% c('PDCD1', 'CD274')] # Remove target gene exp

for(i in 1:5000){ 
  pair <- sample(bc_ISPY_genes, size = 10, replace = FALSE, prob = NULL) 
  null_sr_score <- null_sr.score(ISPY_pemb, pair)
  pval <- null_sr_score[c(which(ISPY_pemb$res_ispy_pemb$pCR == "1"),
                          which(ISPY_pemb$res_ispy_pemb$pCR == "0"))] # Reorder
  flag <- c(rep(1,31),rep(0,38)); target <- target_value(ISPY_pemb)
  target <- target[c(which(ISPY_pemb$res_ispy_pemb$pCR == "1"),
                     which(ISPY_pemb$res_ispy_pemb$pCR == "0"))] # Reorder
  ISPY_bcAUC1[i] <- getAUC(pval, flag)[1] # GI only
  ISPY_bcAUC2[i] <- getAUC((pval* rank.norm(target)), flag)[1] } # BC-SELECT

hist(ISPY_bcAUC1); hist(ISPY_bcAUC2)
round(pnorm(s_ispy_pemb_AUC1, mean = mean(ISPY_bcAUC1), sd = sd(ISPY_bcAUC1), lower.tail=FALSE),4) # GI only
round(pnorm(s_ispy_pemb_AUC2, mean = mean(ISPY_bcAUC2), sd = sd(ISPY_bcAUC2), lower.tail=FALSE),4) # BC-SELECT

# Empirical p-value: target gene expression case
control_bc <-NULL;
control_gene_bc <- ISPY_pemb$genes[-which(ISPY_pemb$genes %in% s_ISPY_pemb_sr_genes_10)] # Remove SR partners 
control_gene_bc <- control_gene_bc[! control_gene_bc %in% c('PDCD1', 'CD274')] # Remove target gene exp

for(i in 1:5000){
  nn_gene <- sample(control_gene_bc, size = 2, replace = FALSE, prob = NULL) 
  pval_bc <- ISPY_pemb$mRNA[nn_gene, c(which(ISPY_pemb$res_ispy_pemb$pCR == "1"), 
                                       which(ISPY_pemb$res_ispy_pemb$pCR == "0"))] # Reorder
  pval_bc2 <- as.numeric(pval_bc[1,])*as.numeric(pval_bc[2,])
  flag <- c(rep(1,31),rep(0,38)); control_bc[i] <- getAUC(rank.norm(pval_bc2), flag)[1]
  rm(nn_gene)}

hist(control_bc) 
round(pnorm(s_ispy_pemb_AUC3, mean = mean(control_bc), sd = sd(control_bc), lower.tail=FALSE),4) # Target


# GSE173839, Durvalumab, ISPY-2
ISPY_bcAUC1 <- NULL; ISPY_bcAUC2 <- NULL # Saved AUC
bc_ISPY_genes <- GSE173839_ISPY2$genes[-which(GSE173839_ISPY2$genes %in% s_ISPY_Durv_sr_genes_10)] # Remove SR partners 
bc_ISPY_genes <- bc_ISPY_genes[! bc_ISPY_genes %in% c('PDCD1', 'CD274')]  # Remove target gene exp

for(i in 1:5000){
  pair <- sample(bc_ISPY_genes, size = 10, replace = FALSE, prob = NULL) 
  null_sr_score <- null_sr.score(GSE173839_ISPY2, pair)
  pval <- null_sr_score[c(which(GSE173839_ISPY2$res_GSE173839_patient$binary == "1"),
                          which(GSE173839_ISPY2$res_GSE173839_patient$binary == "0"))] # Reorder
  flag <- c(rep(1,29),rep(0,42)); target <- target_value(GSE173839_ISPY2)
  target <- target[c(which(GSE173839_ISPY2$res_GSE173839_patient$binary == "1"),
                     which(GSE173839_ISPY2$res_GSE173839_patient$binary == "0"))] # Reorder
  ISPY_bcAUC1[i] <- getAUC(pval, flag)[1] # GI only
  ISPY_bcAUC2[i] <- getAUC( (pval* rank.norm(target)), flag)[1] } # BC-SELECT

hist(ISPY_bcAUC1); hist(ISPY_bcAUC2)
round(pnorm(s_ispy_durv_AUC1, mean = mean(ISPY_bcAUC1), sd = sd(ISPY_bcAUC1), lower.tail=FALSE),4) # GI only
round(pnorm(s_ispy_durv_AUC2, mean = mean(ISPY_bcAUC2), sd = sd(ISPY_bcAUC2), lower.tail=FALSE),4) # BC-SELECT


# Empirical p-value: target gene expression case
control_bc <-NULL;
control_gene_bc <- GSE173839_ISPY2$genes[-which(GSE173839_ISPY2$genes %in% s_ISPY_Durv_sr_genes_10)] # Remove SR partners
control_gene_bc <- control_gene_bc[! control_gene_bc %in% c('PDCD1', 'CD274')] # Remove target gene exp

for(i in 1:5000){
  nn_gene <- sample(control_gene_bc, size = 2, replace = FALSE, prob = NULL) 
  pval_bc <- GSE173839_ISPY2$mRNA[nn_gene, c(which(GSE173839_ISPY2$res_GSE173839_patient$binary == "1"),
                                             which(GSE173839_ISPY2$res_GSE173839_patient$binary == "0"))] # Reorder
  pval_bc2 <- as.numeric(pval_bc[1,])*as.numeric(pval_bc[2,])
  flag <- c(rep(1,29),rep(0,42)); control_bc[i] <- getAUC(rank.norm(pval_bc2), flag)[1]
  rm(nn_gene)}

hist(control_bc)
round(pnorm(s_ispy_durv_AUC3, mean = mean(control_bc), sd = sd(control_bc), lower.tail=FALSE),4) # Target


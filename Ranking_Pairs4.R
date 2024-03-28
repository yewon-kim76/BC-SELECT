# Ranking SL/SDL Pairs
# Load Packages-----------
library(RobustRankAggreg); library(caTools)
library(readxl); library(rmarkdown);library(ggplot2); library(tidyverse); library(dplyr)
library(hrbrthemes); library(viridis); library(forcats); library(data.table); library(Rcpp); 
library(survival); library(parallel); library(pracma); library(ROCR); library(limma);
library(Matrix); library(MKmisc); library(pROC)
#-------------------------

getOR <- function(score,idx.R,idx.NR, thr = 0.44) {
  # function to calculate odds ratio and confidence intervals
  # score is calculated SL scores for all patients
  # idx.R is index of responders in the array of samples
  # idx.NR is index of non responders in the array of samples
  
  rsp <- score[idx.R] # True responders
  non.rsp <- score[idx.NR]  # True non responders
  
  rsp.R <- rsp[rsp > thr] # predicted responders - true positive
  rsp.NR <- rsp[rsp <= thr] # predicted non-responders - false negative
  non.rsp.R <- non.rsp[non.rsp > thr] # predicted responders - false positive
  non.rsp.NR <- non.rsp[non.rsp <= thr] # predicted non-responders - true negative
  
  # Odds ratio
  odds.ratio <- (length(rsp.R)/length(rsp.NR))/(length(non.rsp.R)/length(non.rsp.NR))
  # Standard error of odds ratio
  std.err.or <- sqrt((1/length(rsp.R))+(1/length(rsp.NR))+(1/length(non.rsp.R))+(1/length(non.rsp.NR)))
  
  # Confidence interval - lower bound and upper bound
  lower.bound <- exp(log(odds.ratio) - 1.96 * std.err.or)
  upper.bound <- exp(log(odds.ratio) + 1.96 * std.err.or)
  
  return(c(odds.ratio,lower.bound,upper.bound))
}


#' Rank normalization
#' @param x A numeric \code{vector}
#' @return A normalized \code{vector}
rank_normalization <- function(x) {
  x_back <- x # This just makes sure the NAs are kept in place
  x <- x[!is.na(x)]
  x <- rank(x, ties.method = "average") / length(x)
  x_back[!is.na(x_back)] <- x
  return(x_back) }



#' Bin matrix to values 0, 1 and 2
#' @param x A numeric \code{matrix}
#' @return A binned \code{matrix}
bin_matrix <- function(x, direction = "row", bin = TRUE, split = 3) {
  if (direction == "column") {
    rank_norm <- apply(x, 2, rank_normalization)
    if (split == 2) { q2 <- 1 * (rank_norm >= 1 / 2)}
    if (split == 3) { q2 <- 1 * (rank_norm > 1 / 3) + 1 * (rank_norm > 2 / 3) } }
  if (direction == "row") {
    rank_norm <- t(apply(x, 1, rank_normalization))
    if (split == 2) {q2 <- 1 * (rank_norm >= 1 / 2) }
    if (split == 3) {q2 <- 1 * (rank_norm > 1 / 2) + 1 * (rank_norm > 2 / 3)} }
  if (bin == TRUE) { return(q2)}
  if (bin == FALSE) { return(rank_norm) }  }


getAUC <- function(pval,flag){
        na.inx=which(!is.na(flag) & !is.na(pval))
        pval=pval[na.inx]
        flag=flag[na.inx]
        pred <- prediction(pval, flag)
        perf <- performance(pred,"auc")
        auc=perf@y.values[[1]][1]

        perf1 <- performance(pred, measure="prec", x.measure="rec")
        rec=perf1@x.values[[1]]
        prec=recall=perf1@y.values[[1]]
        prr <- trapz(rec[2:length(rec)], prec[2:length(prec)])
        return(c(auc,prr)) }


# Create a library of candidate partner genes based on single drug targets.
mono_target <- function(drugs, interaction, FDR){
drug_target_map <- readRDS("/gpfs/gsfs12/users/kimy35/SELECT/select_targeted/cox_sort_12_05/BC-SELECT_VERSION CONTROL/drug_target_map.rds") 
gene_targets <- drug_target_map[drug_target_map$drugs == drugs, 2] # Identify drug target gene
  
# Vector of initial candidate gene list
gene_list <- readRDS(file = "/gpfs/gsfs12/users/kimy35/SELECT/select_targeted/cox_sort_12_05/BC-SELECT_VERSION CONTROL/TCGA_genes.rds"); 
gene_list <- gene_list[1:19001] # Remove PD1/PDL1, 19001 genes
  
# Path of processed data
path <- c("/data/kimy35/SELECT/select_targeted/cox_sort_12_05/BC-SELECT_VERSION CONTROL/")
  
# Read essentiality/drug screen/cox analysis/phylogenetic training datasets
ess_file <- grep(list.files(path = paste0(path)), pattern = paste0("^", gene_targets, "_essentiality", 
            collapse = "|"), value = TRUE); ESS <- readRDS(paste0(path, ess_file))
  
drug_file <- grep(list.files(path = paste0(path)), pattern = paste0("^", gene_targets, "_drug_screen", 
             collapse = "|"), value = TRUE); DRUG <- readRDS(paste0(path, drug_file))
  
if(interaction == "SL"){
cox_tcga_file <- grep(list.files(path = paste0(path)), pattern = paste0("^", gene_targets, "_cox_bc", 
                 collapse = "|"), value = TRUE); COX.TCGA <- readRDS(paste0(path, cox_tcga_file))
  
cox_meta_file <- grep(list.files(path = paste0(path)), pattern = paste0("^", gene_targets, "_cox_meta_covariates", 
                 collapse = "|"), value = TRUE); COX.META<- readRDS(paste0(path, cox_meta_file)) }

if(interaction == "SDL"){
cox_tcga_file <- grep(list.files(path = paste0(path)), pattern = paste0("^", gene_targets, "_SDL_cox_bc", 
                   collapse = "|"), value = TRUE); COX.TCGA <- readRDS(paste0(path, cox_tcga_file))
                                                                          
cox_meta_file <- grep(list.files(path = paste0(path)), pattern = paste0("^", gene_targets, "_SDL_cox_meta_covariates", 
                 collapse = "|"), value = TRUE); COX.META<- readRDS(paste0(path, cox_meta_file)) }

phylo_file <- grep(list.files(path = paste0(path)), pattern = paste0("^", gene_targets, "_phylo", collapse = "|"), 
                value = TRUE); PHYLO <- readRDS(paste0(path, phylo_file))
  
# Find SL/SDL candidate pairs
candidates <- data.frame(GENE_A = unlist(lapply(gene_targets, 
                 function(x) rep(x, length(gene_list)))), GENE_B = gene_list)
  
# Essentiality steps
ess_columns <- grep(colnames(ESS), # selecting only relevant columns from essentiality output
                pattern = interaction, # Those with SL, SDL as noted
                value = TRUE, invert = FALSE) # 30 column, ess_split = "tertile", ess_test = "wilcox"
ess_columns <- unique(grep(paste("q2", collapse = "|"), ess_columns, value = TRUE, invert = FALSE))
ess_columns <- unique(grep(paste("linear", collapse = "|"), ess_columns, value = TRUE, invert = TRUE))
  
# Finding lowest p value and adjusted p-value
ess_rank <- ESS %>% # selecting only columns from which to take the union (minimum p-value)
    dplyr::select(all_of(ess_columns)) %>% mutate(union_p_ess = as.numeric(do.call(pmin, c(., na.rm = TRUE)) )) %>%
    as.data.frame() %>% mutate(
    # adding GENE_A and GENE_B information (important for adding to df)
    GENE_A = ESS$GENE_A, GENE_B = ESS$GENE_B) %>%
    subset(GENE_B %in% candidates$GENE_B) %>% # only including candidate pairs in current candidate list
    mutate(fdr_p_ess = p.adjust(.$union_p_ess, method = "BH")) %>%
    # Perform FDR correction (includes all genes, which may also have some NAs, more conservative
    arrange(union_p_ess) %>% # order from lowest to highest p-value only including Gene names and p-values
    dplyr::select(c(GENE_A, GENE_B, union_p_ess, fdr_p_ess))
  
# Drug screen
drug_columns <- grep(colnames(DRUG), # selecting only relevant columns from essentiality output
                pattern = interaction, value = TRUE, invert = FALSE) # Those with SL, SDL as noted, drug_test = "wilcox"
drug_columns <- unique(grep(paste("linear", collapse = "|"), drug_columns,value = TRUE, invert = TRUE)) 
  
# Finding lowest p value and adjusted p-value  
drug_rank <- DRUG %>% # selecting only columns from which to take the union (minimum p-value)
    dplyr::select(all_of(drug_columns)) %>% mutate(union_p_drug = as.numeric(do.call(pmin, c(., na.rm = TRUE)))) %>%
    as.data.frame() %>%  mutate(
    # adding GENE_A and GENE_B information (important for adding to df)
    GENE_A = DRUG$GENE_A, GENE_B = DRUG$GENE_B) %>%
    subset(GENE_B %in% candidates[["GENE_B"]]) %>%  # only including candidate pairs in current candidate list
    mutate(fdr_p_drug = p.adjust(.$union_p_drug, method = "BH")) %>%  # Perform FDR correction
    arrange(union_p_drug) %>% # order from lowest to highest p-value
    dplyr::select(c(GENE_A, GENE_B, union_p_drug, fdr_p_drug))   # only including Gene names and p-values
  
ess_rank <- merge(ess_rank, drug_rank, by = c("GENE_A", "GENE_B"), all = TRUE) # combine data
candidates <- merge(candidates, ess_rank,by = c("GENE_A", "GENE_B"), all.y = TRUE) # sub-setting candidate list
candidates <- candidates %>% subset(fdr_p_ess < FDR | fdr_p_drug < FDR) #fdr_type = "fdr_elm
  
# Cox regression (TCGA + METABRIC): gene_targets
tcga_adj = ifelse(COX.TCGA$coef_adj < 0 , COX.TCGA$pval_adj, 1 ) # preferred negative coefficient
tcga_adj = data.frame(COX.TCGA[,c(2,3,9)],tcga_adj) # adjust covariate and lrt p-values
meta_adj = ifelse(COX.META$coef_adj < 0 , COX.META$pval_adj, 1 ) # preferred negative coefficient
meta_adj = data.frame(COX.META[,c(2,3,9)],meta_adj) # adjust covariate and lrt p-values
cox_result <- merge(tcga_adj, meta_adj, by = c("GENE_A", "GENE_B"), all.x = TRUE, all.y = TRUE)
cox_columns <- colnames(cox_result)[3:6] # 4 related p-values
  
# Finding lowest p-value and adjusted p-value
cox_rank <- cox_result %>% # selecting only columns from which to take the union (minimum p-value)
dplyr::select(all_of(cox_columns)) %>% mutate(union_p_cox = as.numeric(do.call(pmin, c(., na.rm = TRUE)) )) %>%
  as.data.frame() %>% mutate(
  # adding GENE_A and GENE_B information (important for adding to df)
  GENE_A = cox_result$GENE_A, GENE_B = cox_result$GENE_B) %>%
  subset(GENE_B %in% candidates$GENE_B) %>% # only including candidate pairs in current candidate list
  mutate(fdr_p_cox = p.adjust(.$union_p_cox, method = "BH")) %>%
  # Perform FDR correction (includes all genes, which may also have some NAs, more conservative
  arrange(union_p_cox) %>% # order from lowest to highest p-value only including Gene names and p-values
  dplyr::select(c(GENE_A, GENE_B, union_p_cox, fdr_p_cox))
  
cox_rank <- cox_rank[which(cox_rank$fdr_p_cox < FDR),]
candidates <- merge(candidates, cox_rank, by = c("GENE_A", "GENE_B"), all = FALSE)
  
# Phylogenetic
phylo_rank <- PHYLO %>% subset(GENE_B %in% candidates$GENE_B) %>% arrange(phylogenetic_similarity) %>%
  dplyr::select(GENE_A, GENE_B, phylogenetic_similarity) # order from lowest to highest p-value
  
candidates <- merge(candidates, phylo_rank, by = c("GENE_A", "GENE_B"), all = FALSE)
candidates <- candidates %>% subset(phylogenetic_similarity <= 10.5)
  
# Order from lowest to highest p-value of cox analysis
candidates <- candidates[order(candidates$fdr_p_cox),]
return(candidates) }


# Create a library of candidate partner genes based on multiple drug targets.
multiple_target <- function(drugs, interaction, FDR){
drug_target_map <- readRDS("/gpfs/gsfs12/users/kimy35/SELECT/select_targeted/cox_sort_12_05/cox_sort_12_05/BC-SELECT_VERSION CONTROL/drug_target_map.rds") 
gene_targets <- drug_target_map[drug_target_map$drugs == drugs, 2]

# Vector of initial candidate gene list
gene_list <- readRDS(file = "/gpfs/gsfs12/users/kimy35/SELECT/select_targeted/cox_sort_12_05/cox_sort_12_05/BC-SELECT_VERSION CONTROL/TCGA_genes.rds"); 
gene_list <- gene_list[1:19001] # remove PD1/PDL1, 19001 genes

# Path of processed data
path <- c("/data/kimy35/SELECT/select_targeted/cox_sort_12_05/BC-SELECT_VERSION CONTROL/")
RESULT = list();

for(i in 1: length(gene_targets) ){
# Read essentiality/drug screen/cox analysis/phylogenetic training datasets
ess_file <- grep(list.files(path = paste0(path)), pattern = paste0("^", gene_targets[i], "_essentiality", 
            collapse = "|"), value = TRUE); ESS <- readRDS(paste0(path, ess_file))

drug_file <- grep(list.files(path = paste0(path)), pattern = paste0("^", gene_targets[i], "_drug_screen", 
             collapse = "|"), value = TRUE); DRUG <- readRDS(paste0(path, drug_file))

if(interaction == "SL"){
cox_tcga_file <- grep(list.files(path = paste0(path)), pattern = paste0("^", gene_targets[i], "_cox_bc", 
                collapse = "|"), value = TRUE); COX.TCGA <- readRDS(paste0(path, cox_tcga_file))
                                                                          
cox_meta_file <- grep(list.files(path = paste0(path)), pattern = paste0("^", gene_targets[i], "_cox_meta_covariates", 
                collapse = "|"), value = TRUE); COX.META<- readRDS(paste0(path, cox_meta_file)) }

if(interaction == "SDL"){
cox_tcga_file <- grep(list.files(path = paste0(path)), pattern = paste0("^", gene_targets[i], "_SDL_cox_bc", 
                 collapse = "|"), value = TRUE); COX.TCGA <- readRDS(paste0(path, cox_tcga_file))
                                                                          
cox_meta_file <- grep(list.files(path = paste0(path)), pattern = paste0("^", gene_targets[i], "_SDL_cox_meta_covariates", 
                 collapse = "|"), value = TRUE); COX.META<- readRDS(paste0(path, cox_meta_file)) }

phylo_file <- grep(list.files(path = paste0(path)), pattern = paste0("^", gene_targets[i], "_phylo", collapse = "|"), 
                 value = TRUE); PHYLO <- readRDS(paste0(path, phylo_file))

# Find SL/SDL candidate pairs
candidates <- data.frame(GENE_A = unlist(lapply(gene_targets[i], 
                                         function(x) rep(x, length(gene_list)))), GENE_B = gene_list)

# Essentiality steps
ess_columns <- grep(colnames(ESS), # selecting only relevant columns from essentiality output
                    pattern = interaction, # Those with SL, SDL as noted
                    value = TRUE, invert = FALSE) # 30 column, ess_split = "tertile", ess_test = "wilcox"
ess_columns <- unique(grep(paste("q2", collapse = "|"), ess_columns, value = TRUE, invert = FALSE))
ess_columns <- unique(grep(paste("linear", collapse = "|"), ess_columns, value = TRUE, invert = TRUE))

# Finding lowest p value and adjusted p-value
ess_rank <- ESS %>% # selecting only columns from which to take the union (minimum p-value)
  dplyr::select(all_of(ess_columns)) %>% mutate(union_p_ess = as.numeric(do.call(pmin, c(., na.rm = TRUE)) )) %>%
  as.data.frame() %>% mutate(
  # adding GENE_A and GENE_B information (important for adding to df)
  GENE_A = ESS$GENE_A, GENE_B = ESS$GENE_B) %>%
  subset(GENE_B %in% candidates$GENE_B) %>% # only including candidate pairs in current candidate list
  mutate(fdr_p_ess = p.adjust(.$union_p_ess, method = "BH")) %>%
  # Perform FDR correction (includes all genes, which may also have some NAs, more conservative
  arrange(union_p_ess) %>% # order from lowest to highest p-value only including Gene names and p-values
  dplyr::select(c(GENE_A, GENE_B, union_p_ess, fdr_p_ess))

# Drug screen
drug_columns <- grep(colnames(DRUG), # selecting only relevant columns from essentiality output
                pattern = interaction, value = TRUE, invert = FALSE) # Those with SL, SDL as noted, drug_test = "wilcox"
drug_columns <- unique(grep(paste("linear", collapse = "|"), drug_columns, value = TRUE, invert = TRUE)) 

# Finding lowest p-value and adjusted p-value  
drug_rank <- DRUG %>% # selecting only columns from which to take the union (minimum p-value)
  dplyr::select(all_of(drug_columns)) %>% mutate(union_p_drug = as.numeric(do.call(pmin, c(., na.rm = TRUE)))) %>%
  as.data.frame() %>%  mutate(
  # adding GENE_A and GENE_B information (important for adding to df)
  GENE_A = DRUG$GENE_A, GENE_B = DRUG$GENE_B) %>%
  subset(GENE_B %in% candidates[["GENE_B"]]) %>%  # only including candidate pairs in current candidate list
  mutate(fdr_p_drug = p.adjust(.$union_p_drug, method = "BH")) %>%  # Perform FDR correction
  arrange(union_p_drug) %>% # order from lowest to highest p-value
  dplyr::select(c(GENE_A, GENE_B, union_p_drug, fdr_p_drug))   # only including Gene names and p-values

ess_rank <- merge(ess_rank, drug_rank, by = c("GENE_A", "GENE_B"), all = TRUE) # combine data
candidates <- merge(candidates, ess_rank, by = c("GENE_A", "GENE_B"), all.y = TRUE) # subsetting candidate list
candidates <- candidates %>% subset(fdr_p_ess < FDR | fdr_p_drug < FDR) #fdr_type = "fdr_elm

# Cox regression (TCGA + METABRIC): gene_targets
tcga_adj = ifelse(COX.TCGA$coef_adj < 0 , COX.TCGA$pval_adj, 1 ) # preferred negative coefficient
tcga_adj = data.frame(COX.TCGA[,c(2,3,9)],tcga_adj) # adjust covariate and lrt p-values
meta_adj = ifelse(COX.META$coef_adj < 0 , COX.META$pval_adj, 1 ) # preferred negative coefficient
meta_adj = data.frame(COX.META[,c(2,3,9)],meta_adj) # adjust covariate and lrt p-values
cox_result <- merge(tcga_adj, meta_adj, by = c("GENE_A", "GENE_B"), all.x = TRUE, all.y = TRUE)
cox_columns <- colnames(cox_result)[3:6] # 4 related p-values

# Finding lowest p-value and adjusted p-value
cox_rank <- cox_result %>% # selecting only columns from which to take the union (minimum p-value)
  dplyr::select(all_of(cox_columns)) %>% mutate(union_p_cox = as.numeric(do.call(pmin, c(., na.rm = TRUE)) )) %>%
  as.data.frame() %>% mutate(
  # adding GENE_A and GENE_B information (important for adding to df)
  GENE_A = cox_result$GENE_A, GENE_B = cox_result$GENE_B) %>%
  subset(GENE_B %in% candidates$GENE_B) %>% # only including candidate pairs in current candidate list
  mutate(fdr_p_cox = p.adjust(.$union_p_cox, method = "BH")) %>%
  # Perform FDR correction (includes all genes, which may also have some NAs, more conservative
  arrange(union_p_cox) %>% # order from lowest to highest p-value only including Gene names and p-values
  dplyr::select(c(GENE_A, GENE_B, union_p_cox, fdr_p_cox))

cox_rank <- cox_rank[which(cox_rank$fdr_p_cox < FDR),]
candidates <- merge(candidates, cox_rank, by = c("GENE_A", "GENE_B"), all = FALSE)
candidates <- na.omit(candidates)

# Phylogenetic
phylo_rank <- PHYLO %>% subset(GENE_B %in% candidates$GENE_B) %>% arrange(phylogenetic_similarity) %>%
  dplyr::select(GENE_A, GENE_B, phylogenetic_similarity) # order from lowest to highest p-value

candidates <- merge(candidates, phylo_rank, by = c("GENE_A", "GENE_B"), all = FALSE)
candidates <- candidates %>% subset(phylogenetic_similarity <= 10.5)

# Order from lowest to highest p-value of cox analysis
candidates <- candidates[order(candidates$fdr_p_cox),]
RESULT[[i]] <- candidates }

result = bind_rows(lapply(RESULT, bind_rows))
result <- result[order(result$fdr_p_cox),]
return(result) }

# Count up/down regulated genes
down_counts <- function(x){ifelse(x < 0.33, 1, 0)} # check the down-regulated gene for SL
up_counts <- function(x){ifelse(x > (1-0.33), 1, 0)} # check the up-regulated gene for SDL


# AUC calculation with SL genetic interaction
SL_AUC <- function(data1, data2, drugs, drug_response, pair){

drug_target_map <- readRDS("/gpfs/gsfs12/users/kimy35/SELECT/select_targeted/cox_sort_12_05/BC-SELECT_VERSION CONTROL/drug_target_map.rds") 
drug_targets <- drug_target_map[drug_target_map$drugs == drugs, 2]
  
# Calculate SL score based on top nn SL pairs
prediction_data_row_norm <- bin_matrix(data1, direction = "row", bin = FALSE)
prediction_data_col_norm <- bin_matrix(data1, direction = "column", bin = FALSE)
  
sl_score_df <- as.data.frame(matrix(0, nrow = ncol(data1), ncol = 1))
colnames(sl_score_df) <- c("score"); rownames(sl_score_df) <- colnames(prediction_data_row_norm)
potential_sl = active_sl <- 0 # the number of potential pairs

if(sum(!is.na(prediction_data_row_norm[pair, ])) / length(prediction_data_row_norm[pair, ]) > 0.5) {
# IF SL, then if partner is down-regulated, then more likely to respond (increase SL score)
potential_sl <- ifelse(!is.na(prediction_data_row_norm[pair, ]), potential_sl + 1, potential_sl)
active_sl <- ifelse(prediction_data_row_norm[pair, ] < 0.33 & !is.na(prediction_data_row_norm[pair, ]), 
                        active_sl + 1, active_sl) }

# new SL score = (# of inactive pairs / # of total pairs included) in BC360 gene list
prediction_data_row_norm2 <- bin_matrix(data2, direction = "row", bin = FALSE) 
prediction_data_row_norm2 <- na.omit(prediction_data_row_norm2)  # remove NA
whole_down <-  apply(prediction_data_row_norm2, 2 , down_counts)
whole_down_patient <- as.numeric(colSums(whole_down))  # count the down-regulated gene for each patient
whole_sl_score <- whole_down_patient / dim(prediction_data_row_norm2)[1]
  
if(length(drug_targets)  == "1"){ 
  target_gene_exp <- as.numeric(prediction_data_row_norm[drug_targets, ])
  tscore <- (prediction_data_col_norm[drug_targets, ] >= 0.2) * 1 # target gene expression >20% expression
  
  # old SL score = (# of inactive pairs / # of total pairs include) *expression fraction (tscore) of target genes
  sl_score_df$old_sl_score <- (colSums(active_sl, na.rm = TRUE) / colSums(potential_sl,na.rm = TRUE))*tscore 
  
  # new SL score
  sl_score_df$new_sl_score <- (sl_score_df$old_sl_score / whole_sl_score) 
  
  AUC1 <- getAUC(target_gene_exp, drug_response )[1] 
  PR_AUC1 <- getAUC(target_gene_exp, drug_response )[2]
  
  AUC2 <- getAUC(sl_score_df$new_sl_score, drug_response )[1]
  PR_AUC2 <- getAUC(sl_score_df$new_sl_score, drug_response )[2]
                 
  sl_score_df$new_sl_score_target <- (sl_score_df$new_sl_score)*target_gene_exp
  AUC4 <- getAUC(sl_score_df$new_sl_score_target, drug_response )[1] 
  PR_AUC4 <- getAUC(sl_score_df$new_sl_score_target, drug_response )[2]
  
  # empirical p-value: target gene expression case
  control <- NULL
  control_gene <- rownames(data1)[-which(rownames(data1) == drug_targets)] # remove target genes
  
  for(i in 1:5000){
  n_gene <- sample(control_gene, size = 1, replace = FALSE, prob = NULL) 
  control[i] <- getAUC(as.numeric(prediction_data_row_norm[n_gene,]), drug_response )[1]
  rm(n_gene)}
  pval_AUC1 <- pnorm(AUC1, mean = mean(control), sd = sd(control), lower.tail=FALSE) 
  
  # empirical p-value: randomly selected gene case (new SL score)
  bc_metaAUC2 <- NULL; bc_metaAUC4 <- NULL
  bc_meta_gene <- rownames(data1)[-which(rownames(data1) %in% pair)] # remove SL partners 
  bc_meta_gene <- bc_meta_gene[-which(bc_meta_gene == drug_targets)] # remove target genes
  
  for(i in 1:5000){
  npair <- sample(bc_meta_gene, size = nn, replace = FALSE, prob = NULL) 
  nsl_score_df <- as.data.frame(matrix(0, nrow = ncol(data1), ncol = 1))
  colnames(nsl_score_df) <- c("score")
  rownames(nsl_score_df) <- colnames(prediction_data_row_norm)
  potential_sl = active_sl <- 0 # the number of potential pairs
    
  if (sum(!is.na(prediction_data_row_norm[npair, ])) / length(prediction_data_row_norm[npair, ]) > 0.5) {
  # IF SL or SDL_UD, then if partner is down-regulated, then more likely to respond (increase SL score)
  potential_sl <- ifelse(!is.na(prediction_data_row_norm[npair, ]), potential_sl + 1, potential_sl)
  active_sl <- ifelse(prediction_data_row_norm[npair, ] < 0.33 & !is.na(prediction_data_row_norm[npair, ]), 
                      active_sl + 1, active_sl) }
    
  # SL score = # of inactive pairs / # of total pairs include *expression fraction (tscore) of target genes
  nsl_score_df$score <- ((colSums(active_sl,na.rm = TRUE) / colSums(potential_sl,na.rm = TRUE))*tscore) 
  
  bc_metaAUC2[i] <- getAUC(nsl_score_df$score, drug_response )[1] # AUC
  bc_metaAUC4[i] <- getAUC((nsl_score_df$score)*target_gene_exp, drug_response )[1] # AUC
  rm(npair, nsl_score_df )}
  
  pval_AUC2 <- pnorm(AUC2, mean = mean(bc_metaAUC2), sd = sd(bc_metaAUC2), lower.tail=FALSE) 
  pval_AUC4 <- pnorm(AUC4, mean = mean(bc_metaAUC4), sd = sd(bc_metaAUC4), lower.tail=FALSE) } 

if(length(drug_targets) > 1) {
  tscore <- colSums((prediction_data_col_norm[drug_targets, ] >= 0.2) * 1, na.rm = T) / length(drug_targets) 
  # If multiple target genes, then check expression of both
  expression_targets <- prediction_data_row_norm[drug_targets,] # multiple targets
  normal_expression_targets <- as.numeric(colMeans(expression_targets, na.rm = T)) # target expression
  
  # old SL score = (# of inactive pairs / # of total pairs include) *expression fraction (tscore) of target genes
  sl_score_df$old_sl_score <- (colSums(active_sl, na.rm = TRUE) / colSums(potential_sl,na.rm = TRUE))*tscore 
  
  # new SL score
  sl_score_df$new_sl_score <- (sl_score_df$old_sl_score) / whole_sl_score
  
  AUC1 <- getAUC(normal_expression_targets, drug_response )[1]
  PR_AUC1 <- getAUC(normal_expression_targets, drug_response )[2]
  
  AUC2 <- getAUC(sl_score_df$new_sl_score, drug_response )[1]
  PR_AUC2 <- getAUC(sl_score_df$new_sl_score, drug_response )[2]
  
  sl_score_df$new_sl_score_target <- (sl_score_df$new_sl_score)*normal_expression_targets
  AUC4 <- getAUC(sl_score_df$new_sl_score_target, drug_response )[1] 
  PR_AUC4 <- getAUC(sl_score_df$new_sl_score_target, drug_response )[2] 
  
  # empirical p-value: target gene expression case
  control <- NULL
  control_gene <- rownames(data1)[-which(rownames(data1) %in% drug_targets)] # remove target genes
  
  for(i in 1:5000){
  n_gene <- sample(control_gene, size = length(drug_targets), replace = FALSE, prob = NULL) 
  control[i] <- getAUC(as.numeric(colMeans(prediction_data_row_norm[n_gene,])), drug_response)[1]
  rm(n_gene)}
  
  pval_AUC1 <- pnorm(AUC1, mean = mean(control), sd = sd(control), lower.tail=FALSE)  # p-value for target gene expression only
  
  # empirical p-value: randomly selected gene case
  bc_metaAUC2 <- NULL; bc_metaAUC4 <- NULL
  bc_meta_gene <- rownames(data1)[-which(rownames(data1) %in% pair)] # remove SL partners 
  bc_meta_gene <- bc_meta_gene[-which(bc_meta_gene %in% drug_targets)] # remove target genes
  
  for(i in 1:5000){
  npair <- sample(bc_meta_gene, size = nn, replace = FALSE, prob = NULL) 
  nsl_score_df <- as.data.frame(matrix(0, nrow = ncol(data1), ncol = 1))
  colnames(nsl_score_df) <- c("score")
  rownames(nsl_score_df) <- colnames(prediction_data_row_norm)
  potential_sl = active_sl <- 0 # the number of potential pairs
    
  if (sum(!is.na(prediction_data_row_norm[npair, ])) / length(prediction_data_row_norm[npair, ]) > 0.5) {
  # IF SL then if partner is down-regulated, then more likely to respond (increase SL score)
  potential_sl <- ifelse(!is.na(prediction_data_row_norm[npair, ]), potential_sl + 1, potential_sl)
  active_sl <- ifelse(prediction_data_row_norm[npair, ] < 0.33 & !is.na(prediction_data_row_norm[npair, ]), 
                      active_sl + 1, active_sl) }
    
  # SL score = # of inactive pairs / # of total pairs include *expression fraction (tscore) of target genes
  nsl_score_df$score <- ((colSums(active_sl,na.rm = TRUE) / colSums(potential_sl,na.rm = TRUE))*tscore) 
  bc_metaAUC2[i] <- getAUC(nsl_score_df$score, drug_response )[1] ## AUC
  bc_metaAUC4[i] <- getAUC((nsl_score_df$score)*normal_expression_targets, drug_response )[1] ## AUC
  rm(npair, nsl_score_df )}
 
  pval_AUC2 <- pnorm(AUC2, mean = mean(bc_metaAUC2), sd = sd(bc_metaAUC2), lower.tail=FALSE) 
  pval_AUC4 <- pnorm(AUC4, mean = mean(bc_metaAUC4), sd = sd(bc_metaAUC4), lower.tail=FALSE)  } ###

result <- as.data.frame(matrix(c(round(AUC1,3), round(AUC2,3), round(AUC4,3),
                                 round(pval_AUC1,4), round(pval_AUC2,4), round(pval_AUC4,4), 
                                 round(PR_AUC1,3), round(PR_AUC2,3), round(PR_AUC4,3)), nrow = 3, ncol = 3, T))
  colnames(result) <- c("target_exp", "SL","BC-SELECT"); 
  rownames(result) <- c("AUC", "p-val", "PR_AUC")
  return(result)}


# AUC calculation with SDL genetic interaction
SDL_AUC <- function(data1, data2, drugs, drug_response, pair){
  
drug_target_map <- readRDS("/gpfs/gsfs12/users/kimy35/SELECT/select_targeted/cox_sort_12_05/BC-SELECT_VERSION CONTROL/drug_target_map.rds") 
drug_targets <- drug_target_map[drug_target_map$drugs == drugs, 2]
  
# Calculate SDL score based on top nn  pairs
prediction_data_row_norm <- bin_matrix(data1, direction = "row", bin = FALSE)
prediction_data_col_norm <- bin_matrix(data1, direction = "column", bin = FALSE)
  
sdl_score_df <- as.data.frame(matrix(0, nrow = ncol(data1), ncol = 1))
colnames(sdl_score_df) <- c("score"); rownames(sdl_score_df) <- colnames(prediction_data_row_norm)
potential_sdl = active_sdl <- 0 # the number of potential pairs

if(sum(!is.na(prediction_data_row_norm[pair, ])) / length(prediction_data_row_norm[pair, ]) > 0.5) {
# IF SDL_DU, then if partner is up-regulated, then more likely to respond (increase SDL score)
potential_sdl <- ifelse(!is.na(prediction_data_row_norm[pair, ]), potential_sdl + 1, potential_sdl)
active_sdl <- ifelse(prediction_data_row_norm[pair, ] > (1 - 0.33) & !is.na(prediction_data_row_norm[pair, ]), 
                       active_sdl + 1, active_sdl)}

# new SDL score = (# of active pairs / # of total pairs included) in BC360 gene list
prediction_data_row_norm2 <- bin_matrix(data2, direction = "row", bin = FALSE)
prediction_data_row_norm2 <- na.omit(prediction_data_row_norm2)  # remove NA
whole_up <-  apply(prediction_data_row_norm2, 2 , up_counts)
whole_up_patient <- as.numeric(colSums(whole_up))  # count the down-regulated gene for each patient
whole_sdl_score <- whole_up_patient / dim(prediction_data_row_norm2)[1]
  
if(length(drug_targets)  == "1"){ 
  target_gene_exp <- as.numeric(prediction_data_row_norm[drug_targets, ])
  tscore <- (prediction_data_col_norm[drug_targets, ] >= 0.2) * 1 # target gene expression >20% expression
    
  # old SDL score = (# of active pairs / # of total pairs include) *expression fraction (tscore) of target genes
  sdl_score_df$old_sdl_score <- (colSums(active_sdl, na.rm = TRUE) / colSums(potential_sdl, na.rm = TRUE))*tscore 
    
  # new SDL score
  sdl_score_df$new_sdl_score <- (sdl_score_df$old_sdl_score) / whole_sdl_score
    
  AUC1 <- getAUC(target_gene_exp, drug_response )[1] 
  PR_AUC1 <- getAUC(target_gene_exp, drug_response )[2] 
  
  AUC6 <- getAUC(sdl_score_df$new_sdl_score, drug_response )[1]
  PR_AUC6 <- getAUC(sdl_score_df$new_sdl_score, drug_response )[2]
    
  sdl_score_df$new_sdl_score_target <- (sdl_score_df$new_sdl_score)*target_gene_exp
  AUC7 <- getAUC(sdl_score_df$new_sdl_score_target, drug_response )[1] 
  PR_AUC7 <- getAUC(sdl_score_df$new_sdl_score_target, drug_response )[2] 
  
  # empirical p-value: target gene expression case
  control <- NULL
  control_gene <- rownames(data1)[-which(rownames(data1) == drug_targets)] # remove target genes
    
  for(i in 1:5000){
    n_gene <- sample(control_gene, size = 1, replace = FALSE, prob = NULL) ##
    control[i] <- getAUC(as.numeric(prediction_data_row_norm[n_gene,]), drug_response )[1]
    rm(n_gene)}
  pval_AUC1 <- pnorm(AUC1, mean = mean(control), sd = sd(control), lower.tail=FALSE) # p-value for target gene expression only
    
  # empirical p-value: randomly selected gene case
  bc_metaAUC6 <- NULL; bc_metaAUC7 <- NULL
  bc_meta_gene <- rownames(data1)[-which(rownames(data1) %in% pair)] # remove SDL partners 
  bc_meta_gene <- bc_meta_gene[-which(bc_meta_gene == drug_targets)] # remove target genes
    
  for(i in 1:5000){
  npair <- sample(bc_meta_gene, size = nn, replace = FALSE, prob = NULL) 
  nsdl_score_df <- as.data.frame(matrix(0, nrow = ncol(data1), ncol = 1))
  colnames(nsdl_score_df) <- c("score")
  rownames(nsdl_score_df) <- colnames(prediction_data_row_norm)
  potential_sdl = active_sdl <- 0 # the number of potential pairs
      
  if (sum(!is.na(prediction_data_row_norm[npair, ])) / length(prediction_data_row_norm[npair, ]) > 0.5) {
  # IF SDL_UD, then if partner is up-regulated, then more likely to respond (increase SDL score)
  potential_sdl <- ifelse(!is.na(prediction_data_row_norm[npair, ]), potential_sdl + 1, potential_sdl)
  active_sdl <- ifelse(prediction_data_row_norm[npair, ] > (1 - 0.33) & !is.na(prediction_data_row_norm[npair, ]), 
                          active_sdl + 1, active_sdl) }
      
  # SDL score = # of active pairs / # of total pairs include *expression fraction (tscore) of target genes
  nsdl_score_df$score <- ((colSums(active_sdl,na.rm = TRUE) / colSums(potential_sdl,na.rm = TRUE))*tscore) 
  bc_metaAUC6[i] <- getAUC(nsdl_score_df$score, drug_response )[1] # AUC
  bc_metaAUC7[i] <- getAUC((nsdl_score_df$score)*target_gene_exp, drug_response )[1] # AUC
  rm(npair, nsdl_score_df )}

  pval_AUC6 <- pnorm(AUC6, mean = mean(bc_metaAUC6), sd = sd(bc_metaAUC6), lower.tail=FALSE) 
  pval_AUC7 <- pnorm(AUC7, mean = mean(bc_metaAUC7), sd = sd(bc_metaAUC7), lower.tail=FALSE) }
  
  
  if(length(drug_targets) > 1) {
  tscore <- colSums((prediction_data_col_norm[drug_targets, ] >= 0.2) * 1, na.rm = T) / length(drug_targets) 
  # If multiple target genes, then check expression of both
  expression_targets <- prediction_data_row_norm[drug_targets,] # multiple targets
  normal_expression_targets <- as.numeric(colMeans(expression_targets, na.rm = T)) # target expression
    
  # old SDL score = (# of active pairs / # of total pairs include) *expression fraction (tscore) of target genes
  sdl_score_df$old_sdl_score <- (colSums(active_sdl, na.rm = TRUE) / colSums(potential_sdl,na.rm = TRUE))*tscore 
    
  # new SDL score
  sdl_score_df$new_sdl_score <- (sdl_score_df$old_sdl_score) / whole_sdl_score
    
  AUC1 <- getAUC(normal_expression_targets, drug_response)[1]
  PR_AUC1 <- getAUC(normal_expression_targets, drug_response)[2]
  
  AUC6 <- getAUC(sdl_score_df$new_sdl_score, drug_response)[1]
  PR_AUC6 <- getAUC(sdl_score_df$new_sdl_score, drug_response)[2]
    
  sdl_score_df$new_sdl_score_target <- (sdl_score_df$new_sdl_score)*normal_expression_targets
  AUC7 <- getAUC(sdl_score_df$new_sdl_score_target, drug_response)[1] 
  PR_AUC7 <- getAUC(sdl_score_df$new_sdl_score_target, drug_response)[2] 
  
  # empirical p-value: target gene expression case
  control <- NULL
  control_gene <- rownames(data1)[-which(rownames(data1) %in% drug_targets)] # remove target genes
    
  for(i in 1:5000){
  n_gene <- sample(control_gene, size = length(drug_targets), replace = FALSE, prob = NULL) 
  control[i] <- getAUC(as.numeric(colMeans(prediction_data_row_norm[n_gene,])), drug_response)[1]
  rm(n_gene)}
    
  pval_AUC1 <- pnorm(AUC1, mean = mean(control), sd = sd(control), lower.tail=FALSE) # p-value for target gene expression only
    
  # empirical p-value: randomly selected gene case
  bc_metaAUC6 <- NULL; bc_metaAUC7 <- NULL
  bc_meta_gene <- rownames(data1)[-which(rownames(data1) %in% pair)] # remove SDL partners 
  bc_meta_gene <- bc_meta_gene[-which(bc_meta_gene %in% drug_targets)] # remove target genes
    
  for(i in 1:5000){
  npair <- sample(bc_meta_gene, size = nn, replace = FALSE, prob = NULL) 
  nsdl_score_df <- as.data.frame(matrix(0, nrow = ncol(data1), ncol = 1))
  colnames(nsdl_score_df) <- c("score")
  rownames(nsdl_score_df) <- colnames(prediction_data_row_norm)
  potential_sdl = active_sdl <- 0 # the number of potential pairs
      
  if (sum(!is.na(prediction_data_row_norm[npair, ])) / length(prediction_data_row_norm[npair, ]) > 0.5) {
  # IF SDL_DU, then if partner is up-regulated, then more likely to respond (increase SDL score)
  potential_sdl <- ifelse(!is.na(prediction_data_row_norm[npair, ]), potential_sdl + 1, potential_sdl)
  active_sdl <- ifelse(prediction_data_row_norm[npair, ] > (1 - 0.33) & !is.na(prediction_data_row_norm[npair, ]), 
                            active_sdl + 1, active_sdl) }
      
  # SDL score = # of active pairs / # of total pairs include *expression fraction (tscore) of target genes
  nsdl_score_df$score <- ((colSums(active_sdl,na.rm = TRUE) / colSums(potential_sdl,na.rm = TRUE))*tscore) 
  bc_metaAUC6[i] <- getAUC(nsdl_score_df$score, drug_response )[1] ## AUC
  bc_metaAUC7[i] <- getAUC((nsdl_score_df$score)*normal_expression_targets, drug_response )[1] ## AUC
  rm(npair, nsdl_score_df )}

  pval_AUC6 <- pnorm(AUC6, mean = mean(bc_metaAUC6), sd = sd(bc_metaAUC6), lower.tail=FALSE) 
  pval_AUC7 <- pnorm(AUC7, mean = mean(bc_metaAUC7), sd = sd(bc_metaAUC7), lower.tail=FALSE) }
  
  result <- as.data.frame(matrix(c(round(AUC1,3), round(AUC6,3), round(AUC7,3), 
                                   round(pval_AUC1,4), round(pval_AUC6,4),round(pval_AUC7,4), 
                                   round(PR_AUC1,3), round(PR_AUC6,3), round(PR_AUC7,3)), nrow = 3, ncol = 3, T))
  colnames(result) <- c("target_exp", "SDL","BC-SELECT"); 
  rownames(result) <- c("AUC", "p-val", "PR_AUC")
  
  return(result) }

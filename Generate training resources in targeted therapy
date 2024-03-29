rm(list=ls())
# Load libraries -----------------------------------------------------
library(optparse);library(parallel);library(tidyverse);library(survival)
library(Rcpp)

# Essentiality
# Related functions ----------------------------------------------
molecdata_midpoint <- function(molecdata) {
  median_molecdata <- median(molecdata, na.rm = TRUE)
  if ((median_molecdata == 1 | median_molecdata == 0 | is.na(median_molecdata))) {
    median_molecdata <- 0.5}
  median_molecdata }

split_at_midpoint <- function(molecdata, fitness) {
  median_molecdata <- molecdata_midpoint(molecdata)
  low <- fitness[molecdata <= median_molecdata]
  high <- fitness[molecdata > median_molecdata]
  list(low = low, high = high)}

split_one_two <- function(fitness, molecdata_q2) {
  zero <- fitness[molecdata_q2 == 0]
  two <- fitness[molecdata_q2 == 2]
  list(zero = zero, two = two)}

bidirectional_wilcox_test <- function(x, y) {
  less <- wilcox_test_na(x, y, alternative1 = "less")
  greater <- wilcox_test_na(x, y, alternative1 = "greater")
  list(less = less, greater = greater)}

lm_summary_coef <- function(form = forumlaString, data) {
  tryCatch(
    {model <- summary(lm(as.formula(form), data))
      if (any(model$aliased)) {
        matrix(as.numeric(NA),
               nrow = length(model$aliased),
               ncol = ncol(model$coefficients),
               dimnames = list(names(model$aliased), colnames(model$coefficients)) )
      } else {
        model$coefficients}},
    error = function(e) {
      matrix(NA, nrow = 4, ncol = 4)} ) }

wilcox_test_na <- function(x, y, alternative1, paired = FALSE) { # Wilcox test, returns p-value
  tryCatch(
    wilcox.test( x, y, alternative = alternative1, paired = paired)$p.value,
    error = function(e) NA ) }

# Functions --------------------------------------------------------------------
essentiality <- function(screen = screen_data, # R list object with several data slots
                         sl_pair = sl_pair, 
                         data_slots = c("mrna"), # Remove scna (copynumber removed)
                         datasets = "", # Utilize all objects
                         q2 = FALSE) {
  
  # Results columns
  results <- list()

  # Just in case make sure the sl_pair does not contain NA values
  if (sum(is.na(sl_pair)) == 0) {
    gene_a <- sl_pair["GENE_A"]; gene_b <- sl_pair["GENE_B"]
    
  results["GENE_A"] <- gene_a; results["GENE_B"] <- gene_b
  for (set in datasets) {
    fitness <- screen[[set]][["gene_fitness"]] # for each data modality
    for (slot in data_slots) {
    # Get essentiality/gene fitness (sometimes also referred to gene effect score for A -> B
    model_data <- screen[[set]][["meta"]]
    molcdata <- screen[[set]][[slot]]
        
    model_data$GENE_A_fit <- fitness[gene_a, ]
    model_data$GENE_A_exp <- molcdata[gene_a, ]
    model_data$GENE_B_fit <- fitness[gene_b, ]
    model_data$GENE_B_exp <- molcdata[gene_b, ]
        
    # Unadjusted A -> B
    formulaString1 <- paste("GENE_B_fit ~ GENE_A_exp")
    a_b_lm_unadj <- lm_summary_coef(form = formulaString1, data = model_data)
    coef1_unadj <- a_b_lm_unadj[2, 1]; pval1_unadj <- a_b_lm_unadj[2, 4]
        
    results[paste("SL_A_B_SR_UD_linear_unadj", slot, set, sep = "_")] <-
    ifelse(coef1_unadj > 0, pval1_unadj, ifelse(is.na(coef1_unadj), NA, 1))
        
    results[paste("SDL_UD_SR_DD_linear_unadj", slot, set, sep = "_")] <-
    ifelse(coef1_unadj < 0, pval1_unadj, ifelse(is.na(coef1_unadj), NA, 1))
    
    formulaString3 <- paste("GENE_A_fit ~ GENE_B_exp")
        b_a_lm_unadj <- lm_summary_coef( form = formulaString3, data = model_data)
        coef2_unadj <- b_a_lm_unadj[2, 1]; pval2_unadj <- b_a_lm_unadj[2, 4]
        
    results[paste("SL_B_A_SR_DU_linear_unadj", slot, set, sep = "_")] <-
          ifelse(coef2_unadj > 0, pval2_unadj, ifelse(is.na(coef2_unadj), NA, 1))
    
    results[paste("SDL_DU_SR_DD_linear_unadj", slot, set, sep = "_")] <-
          ifelse(coef2_unadj < 0, pval2_unadj, ifelse(is.na(coef2_unadj), NA, 1))
    
    # This is imporatnt if looking at a binarized fusion gene or something where the "expression" is a 0 or 1
    
    results[[paste("GENE_A", set, "fit_score#cl", sep = "_")]] <- sum(!is.na(fitness[gene_a, ]))
    results[[paste("GENE_B", set, "fit_score#cl", sep = "_")]] <- sum(!is.na(fitness[gene_b, ]))
        
    fit_ab <- split_at_midpoint(
          molecdata = screen[[set]][[slot]][gene_a, ], fitness = fitness[gene_b, ])
        
    test <- bidirectional_wilcox_test(fit_ab$low, fit_ab$high)
    results[paste("SL_A_B_SR_UD_wilcox", slot, set, sep = "_")] <- test$less
    results[paste("SDL_UD_SR_DD_wilcox", slot, set, sep = "_")] <- test$greater
        
    # Change direction B -> A
    fit_ba <- split_at_midpoint(
          molecdata = screen[[set]][[slot]][gene_b, ], fitness = fitness[gene_a, ])
        
    test <- bidirectional_wilcox_test(fit_ba$low, fit_ba$high)
    results[paste("SL_B_A_SR_DU_wilcox", slot, set, sep = "_")] <- test$less
    results[paste("SDL_DU_SR_DD_wilcox", slot, set, sep = "_")] <- test$greater} 
    
    if (q2) { q2_slots <- paste0(data_slots, "_q2") # q2 slots
       for (slot in q2_slots) {
        b_given_a <- split_one_two(
          fitness = fitness[gene_b, ], molecdata_q2 = screen[[set]][[slot]][gene_a, ])
        
        test <- bidirectional_wilcox_test(b_given_a$zero, b_given_a$two)
        results[paste("SL_A_B_SR_UD_wilcox", slot, set, sep = "_")] <- test$less
        results[paste("SDL_UD_SR_DD_wilcox", slot, set, sep = "_")] <- test$greater
        
        # Change direction B -> A
        a_given_b <- split_one_two(
          fitness = fitness[gene_a, ], molecdata_q2 = screen[[set]][[slot]][gene_b, ])
        
        test <- bidirectional_wilcox_test(a_given_b$zero, a_given_b$two)
        results[paste("SL_B_A_SR_DU_wilcox", slot, set, sep = "_")] <- test$less
        results[paste("SDL_DU_SR_DD_wilcox", slot, set, sep = "_")] <- test$greater }}
    
    } }
   return(results)}

# fit--------------------------------------------------------------------
screen <- readRDS("/gpfs/gsfs12/users/kimy35/SELECT/select_targeted/data/new_essentiality_object_joo.rds")
drug_target = "CYP19A1" # Define drug targets
ess_result <- NULL;

for(i in 1:length(screen$prob_ach_new2$genes)){
  gene <- screen$prob_ach_new2$genes[i] # All candidates of partner genes
  test <- essentiality(screen = screen,
                      sl_pair = c(GENE_A = drug_target, GENE_B = gene),
                      data_slots = c("mrna"), # Remove copynumber
                      datasets = c("prob_ach_new2", "prob_ach_new1", "prob_ach_old1", 
                                   "prob_mar_new1", "prob_mar_old1"),
                      q2 = TRUE)
  ess_result <- rbind(ess_result, data.frame(test))
  print(i) }

setwd("/gpfs/gsfs12/users/kimy35/SELECT/select_targeted/Target_training_BC-SELECT")
saveRDS(ess_result, file = "CYP19A1_essentiality.rds") # Save the results

rm(gene, i, test)

# Drug target (continued from Essentiality)
# Related functions ----------------------------------------------
min_na <- function(values) {
  ifelse(all(is.na(values)), as.numeric(NA), min(values, na.rm = TRUE)) }

# Functions --------------------------------------------------------------------
drug_screen <- function(screen = screen_data,
                        sl_pair = sl_pairs,
                        data_slots = "mrna",
                        datasets = "") {
  
  # results columns
  results <- list()
  
  # Just in case make sure the sl_pair does not contain NA values
  if (sum(is.na(sl_pair)) == 0) {
    results[["GENE_A"]] <- sl_pair["GENE_A"]
    results[["GENE_B"]] <- sl_pair["GENE_B"]
    
  for (set in datasets) {
    # For each data modality
    drugs <- screen[[set]][["drug_targets"]][, "drugs"]
    gene_targets <- screen[[set]][["drug_targets"]][, "gene_targets"]
    drugs_target_gene_a <- drugs[gene_targets == sl_pair[["GENE_A"]]]
    drugs_target_gene_b <- drugs[gene_targets == sl_pair[["GENE_B"]]]
      
    for (slot in data_slots) {
      slot_q2_matrix <- screen[[set]][[paste0(slot, "_q2")]]
      # Use the q2 row for gene b to split for wilcox A -> B.
      molecdata <- slot_q2_matrix[sl_pair["GENE_B"], ]
        
      drug_df_sl_a_b <- double(); drug_df_sdl_du <- double()
      for (drug in drugs_target_gene_a) {
        # Get essentiality/gene fitness (sometimes also referred to gene effect score for A -> B)
        fitness <- screen[[set]][["gene_fitness"]][drug, , drop = TRUE]
        pval <- bidirectional_wilcox_test(fitness[molecdata <= 0], fitness[molecdata >= 1])
          
        drug_df_sl_a_b[paste("SL_A_B_SR_DU", slot, set, drug, sep = "_")] <- pval$less
        drug_df_sdl_du[paste("SDL_DU_SR_DD", slot, set, drug, sep = "_")] <- pval$greater }
        
      # Split for wilcox B -> A
      molecdata <- slot_q2_matrix[sl_pair["GENE_A"], ]
        
      drug_df_sl_b_a <- double(); drug_df_sdl_ud <- double()
      for (drug in drugs_target_gene_b) {
        fitness <- screen[[set]][["gene_fitness"]][drug, , drop = TRUE]
        pval <- bidirectional_wilcox_test(fitness[molecdata <= 0], fitness[molecdata >= 1])
          
        drug_df_sl_b_a[paste("SL_B_A_SR_DU", slot, set, drug, sep = "_")] <- pval$less
        drug_df_sdl_ud[paste("SDL_UD_SR_DD", slot, set, drug, sep = "_")] <- pval$greater }
        
        results[[paste("SL_A_B_SR_DU", slot, set, sep = "_")]] <- min_na(drug_df_sl_a_b)
        results[[paste("SDL_DU_SR_DD", slot, set, sep = "_")]] <- min_na(drug_df_sdl_du)
        results[[paste("SL_B_A_SR_DU", slot, set, sep = "_")]] <- min_na(drug_df_sl_b_a)
        results[[paste("SDL_UD_SR_DD", slot, set, sep = "_")]] <- min_na(drug_df_sdl_ud) } } }
  return(results) }


# fit--------------------------------------------------------------------
drug_screen_result <- NULL;
drug_screen_data <- readRDS("/gpfs/gsfs12/users/kimy35/SELECT/select_targeted/data/new_drug_object_joo.rds")

for(i in 1:length(drug_screen_data$basu$genes)){
  gene <- drug_screen_data$basu$genes[i] # All candidates of partner genes
  test <- drug_screen( screen = drug_screen_data,
                       sl_pair = c(GENE_A = drug_target, GENE_B = gene),
                       data_slots = c("mrna"), 
                       datasets = c("basu", "gdsc", "ccle") )
  drug_screen_result <- rbind(drug_screen_result, data.frame(test))
  print(i) }

setwd("/gpfs/gsfs12/users/kimy35/SELECT/select_targeted/Target_training_BC-SELECT")
saveRDS(drug_screen_result, file = "CYP19A1_drug_screen.rds") # Save the results

rm(i, gene, test, screen, drug_screen_data)


# Cox regression
# Related functions ----------------------------------------------
qnorm_normalization <- function(x) {
  # Saving x_back allows us to keep NAs that rank() removes.
  x_back <- x
  x <- rank(x, ties.method = "average", na.last = NA)
  x <- qnorm(x / (length(x) + 1))
  x_back[!is.na(x_back)] <- x
  return(x_back)}

# Functions -------------------------------------------------------------------
cox_sl_pair <- function(sl_pair = c(GENE_A = "ALK", GENE_B = "LINC00115"), # Define genes
                        prob = data, # Defind RNA-seq data
                        clinical = clinical, # Defind phenotype data
                        interaction = "SL",
                        survival = "survival"){
  rows <- list()
  
  if (sl_pair["GENE_A"] %in% prob$genes & sl_pair["GENE_B"] %in% prob$genes) {
  # Only include cases where SL pairs are in dataset
  # Build dataframe for cox regression (all variables and covariates)
  cox_df <- clinical %>%
    add_column(
      g1 = prob$mrna_qnorm[sl_pair["GENE_A"], ], g2 = prob$mrna_qnorm[sl_pair["GENE_B"], ],
      g1_q2 = prob$mrna_q2[sl_pair["GENE_A"], ], g2_q2 = prob$mrna_q2[sl_pair["GENE_B"], ] )
  
    cox_df <- cox_df[complete.cases(cox_df), ] # Only include samples where all covariates and cov +g1 +g2 are non NA

  if (nrow(cox_df) > 10) { # Must have at least 10 patients with all variables going into model
     
    cox_df <- cox_df %>% mutate(cov = as.numeric(g1_q2 == 0 & g2_q2 == 0)) # Both of them are lower expressed, cov = 1
      
    # Model 1: Model without controlling for confounding factors
    uncntrl <- coxph(Surv(time, status) ~ cov, data = cox_df)
    model1 <- summary(uncntrl)

    # Add possible covariates
    var <- cox_df %>% select(any_of(c("types","age","sex","race","stage")))
    f1.list = names(var) # Select any possible variables
    f1.list <- ifelse(f1.list %in% c("types","sex","race","stage"), paste("strata(",f1.list,")", sep="" ), f1.list) # change to strata
      
    # Model 2: Reduced model, without the interaction term (cov)
    formulaString <- paste( "Surv(time, status) ~", paste(f1.list, collapse = " + ") )
    cntrl_reduced <- coxph(as.formula(formulaString), data = cox_df)
    model2 <- summary(cntrl_reduced)
     
    # Model 3: Full model, with interaction term (cov)
    formulaString <- paste( "Surv(time, status) ~ cov +", paste(f1.list, collapse = " + ") )
    cntrl_full <- coxph(as.formula(formulaString), data = cox_df)
    model3 <- summary(cntrl_full)
   
 # In the following tests df, Save coefficients and p-values
 # lrt = Log likelihood ratio test. To test the significance of the interaction 
 # We need to compute the difference between the log likelihood statistic of the reduced model (model 2) and full model (model 3)

    cov_lrt <- (2) * (model3$loglik[2] - model2$loglik[2])
    res <- list(
      test = "cov",
      GENE_A = drug_target,
      GENE_B = sl_pair["GENE_B"],
      coef_unadj = model1$coefficients[1], # cov coefficient, without covariates model
      pval_unadj = model1$coefficients[1, 5], # cov p-value, without covariates model
      coef_adj = model3$coefficients[1], # cov coefficient, with covariates model
      pval_adj = model3$coefficients[1, 5], # cov p-value, with covariates model
      lrt = cov_lrt,
      pval_lrt = 1 - pchisq(cov_lrt, df = 1),
      pts_w = sum(cox_df$cov == 1),
      pts_wo = sum(cox_df$cov == 0),
      events_w = sum(cox_df$status == 1 & cox_df$cov == 1),
      events_wo = sum(cox_df$status == 1 & cox_df$cov == 0))
      
    rows[[length(rows) + 1]] <- res } }
  
  bind_rows(rows) }

# fit--------------------------------------------------------------------
# Patient level data: TCGA and METABRIC
load("/data/kimy35/prob.curtis.extended.RData") # METABRIC: metabric.list, metabric.surv
load("/gpfs/gsfs12/users/kimy35/SELECT/select_immuno/data/prob.TCGA.ICI.RData") # TCGA breast cancer

cox_bc_results <- NULL; 
for(i in 1:19001){
  gene = TCGA.list$genes[i] # All candidates of partner genes
  test_TCGA <- cox_sl_pair(sl_pair = c(GENE_A = drug_target, GENE_B = gene),
                           prob = BC.list, # Defind RNA-seq data
                           clinical = BC.surv, # Defind phenotype data
                           interaction = "SL",
                           survival = "survival")
  cox_bc_results <- rbind(cox_bc_results, test_TCGA)
  print(i)}


setwd("/gpfs/gsfs12/users/kimy35/SELECT/select_targeted/Target_training_BC-SELECT")
saveRDS(cox_bc_results, file = "CYP19A1_cox_bc.rds") # Save the results

rm(prob, clinical , cntrl_full, cntrl_reduced, cox_df, model1, model2, model3, res, uncntrl, var, 
   cov_lrt, f1.list, formulaString, lrt, pval_lrt, gene, i, test_TCGA)


# Phylogenetic Matching
# Functions -------------------------------------------------------------------
# Identify phylogenetically linked candidates of SL pair (high phylogenetic similarity) by analyzing phylogenetic profiles.
# In a candidate SL pair across 86 species (adopting the method of Tabach et al.40,41).
# We then calculate the phylogenetic similarity between gene A and gene B using a non-negative matrix
# using a non-negative matrix factorization (NMF).
# phylogenetic similarity between A and B using a non-negative matrix, which measures the Euclidian distance while
# taking into account their phylogeny. To determine the threshold used for this step,
# we used the median phylogenetic similarity of the pairs in Initial Set I
# as it optimally separates the positive and negative sets among these pairs.

# Phylogenetic profiles of 19017 human genes across 86 eukaryotic genomes.
# The matrix was normalized by the evolutionary distances between organisms
# and the protein length and clustered using average linkage.
# The entry values are between 0 and 1 where 1 represents
# 100% identity and 0 corresponds to no detectable homolog.
# The phylogenetic profile is downloaded from Yuval Tabach et al., Mol Syst Biol. (2013), Supplementary Table 1

# feature_weight Matrix: The feature weights are determined based on the phylogenetic tree 
#    (Ensembl database: http://useast.ensembl.org/index.html)

phylogenetic_profiles <- function(phylo_profiles = phylo_profiles,
                                  feature_weight = feature_weight,
                                  sl_pairs) {

  sl_phylo <- cbind(  # Get ids
    match(sl_pairs[1], phylo_profiles$genes),
    match(sl_pairs[2], phylo_profiles$genes))
  
  # Get phylo profile for gene A and gene B
  phylo_geneA <- phylo_profiles[sl_phylo[, 1], -(1:3)]
  phylo_geneB <- phylo_profiles[sl_phylo[, 2], -(1:3)]
  
  # Get the phylogenetic distance between gene A and B
  featureMat <- (phylo_geneA - phylo_geneB)^2
  featureMat <- data.matrix(featureMat)
  
  # %*% is equivalent to multiplying to matrices then sum up the products
  nmf_res <- featureMat %*% t(feature_weight)
  nmf_res <- data.frame( t(sl_pairs), nmf_res[1], stringsAsFactors=FALSE)
  colnames(nmf_res) <- c("GENE_A", "GENE_B", "phylogenetic_similarity")
  nmf_res$rank_norm <- rank(
    nmf_res$phylogenetic_similarity,
    na.last = "keep" ) / length(nmf_res$phylogenetic_similarity)
  
  return(nmf_res) }

# fit--------------------------------------------------------------------
phylo_profiles <- readRDS("/gpfs/gsfs12/users/kimy35/SELECT/select_targeted/data/phylo.rds")
feature_weight <- readRDS("/gpfs/gsfs12/users/kimy35/SELECT/select_targeted/data/feature_weight.rds")

phylo_result <- NULL 

for(i in 1:length(phylo_profiles$genes)){
  gene = as.character(phylo_profiles$genes)[i]
  sl_pairs = c(drug_target, gene)
  
  phylogenetic_similarity <- phylogenetic_profiles(phylo = phylo_profiles,
                                                   feature_weight = feature_weight,
                                                   sl_pairs = sl_pairs)
  phylo_result <-rbind(phylo_result, phylogenetic_similarity )
  print(i)}

setwd("/gpfs/gsfs12/users/kimy35/SELECT/select_targeted/Target_training_BC-SELECT")
saveRDS(phylo_result, file = "CYP19A1_phylo.rds") # Save the results

rm(i,gene, sl_pairs, phylogenetic_similarity)

# How to select? (Lee et al. Nat Comm)
# Select candidate SL pairs that pass have small phylogenetic distance measured by NMF (< 10.5).
# ISLE uses this criterion to determine phylogenetically linked pairs.
# This threshold will be applied to identify SL partner genes.

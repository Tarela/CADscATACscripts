#Load packages
library(tidyverse)
library(data.table)

#Read in rasqual results
results <- fread("Endothelial_rasqual_results_allpeaks_022421.txt")
results <- as_tibble(results)

#Read in permuted rasqual results
results_permuted1 <- fread("Endothelial_rasqual_results_allpeaks_perm1.txt")
results_permuted1 <- as_tibble(results_permuted1)
results_permuted2 <- fread("Endothelial_rasqual_results_allpeaks_perm2.txt")
results_permuted2 <- as_tibble(results_permuted2)
results_permuted3 <- fread("Endothelial_rasqual_results_allpeaks_perm3.txt")
results_permuted3 <- as_tibble(results_permuted3)
results_permuted4 <- fread("Endothelial_rasqual_results_allpeaks_perm4.txt")
results_permuted4 <- as_tibble(results_permuted4)
results_permuted5 <- fread("Endothelial_rasqual_results_allpeaks_perm5.txt")
results_permuted5 <- as_tibble(results_permuted5)


#Make column names
colnames(results) <- c("Feature_ID", "rsID", "Chromosome", "SNP_position", "Ref_allele",
                       "Alt_allele", "Freq", "HWE_Chisquare", "Imp_quality", "Log10_BH_Q",
                       "Chisquare", "Effect_size", "Delta", "Phi", "Overdispersion",
                       "SNP_id_region", "Num_feature_SNPs", "Num_tested_SNPs", "Num_iterations_null", "Num_iterations_alt",
                       "Random_ties", "Log_likelihood_null", "Convergence_status", "Sq_corr_fSNPs", "Sq_corr_rSNP")

colnames(results_permuted1) <- colnames(results)
colnames(results_permuted2) <- colnames(results)
colnames(results_permuted3) <- colnames(results)
colnames(results_permuted4) <- colnames(results)
colnames(results_permuted5) <- colnames(results)


#Remove skipped rows
results <- results %>% filter(rsID != "SKIPPED")
results_permuted1 <- results_permuted1 %>% filter(rsID != "SKIPPED")
results_permuted2 <- results_permuted2 %>% filter(rsID != "SKIPPED")
results_permuted3 <- results_permuted3 %>% filter(rsID != "SKIPPED")
results_permuted4 <- results_permuted4 %>% filter(rsID != "SKIPPED")
results_permuted5 <- results_permuted5 %>% filter(rsID != "SKIPPED")


#Calculate p-values
results <- results %>% mutate(P = pchisq(Chisquare, 1, lower = F))
results_permuted1 <- results_permuted1 %>% mutate(P = pchisq(Chisquare, 1, lower = F))
results_permuted2 <- results_permuted2 %>% mutate(P = pchisq(Chisquare, 1, lower = F))
results_permuted3 <- results_permuted3 %>% mutate(P = pchisq(Chisquare, 1, lower = F))
results_permuted4 <- results_permuted4 %>% mutate(P = pchisq(Chisquare, 1, lower = F))
results_permuted5 <- results_permuted5 %>% mutate(P = pchisq(Chisquare, 1, lower = F))

#Sort by p-value
results_sorted <- arrange(results, P)
results_permuted1_sorted <- arrange(results_permuted1, P)
results_permuted2_sorted <- arrange(results_permuted2, P)
results_permuted3_sorted <- arrange(results_permuted3, P)
results_permuted4_sorted <- arrange(results_permuted4, P)
results_permuted5_sorted <- arrange(results_permuted5, P)

#Convert Log10q to q-value
results_sorted <- results_sorted %>% mutate(q = 10^(Log10_BH_Q))
results_permuted1_sorted <- results_permuted1_sorted %>% mutate(q = 10^(Log10_BH_Q))
results_permuted2_sorted <- results_permuted2_sorted %>% mutate(q = 10^(Log10_BH_Q))
results_permuted3_sorted <- results_permuted3_sorted %>% mutate(q = 10^(Log10_BH_Q))
results_permuted4_sorted <- results_permuted4_sorted %>% mutate(q = 10^(Log10_BH_Q))
results_permuted5_sorted <- results_permuted5_sorted %>% mutate(q = 10^(Log10_BH_Q))

#Extract q values from each permutation run
perm_df <- data.frame(results_permuted1_sorted$q, results_permuted2_sorted$q, results_permuted3_sorted$q, results_permuted4_sorted$q, results_permuted5_sorted$q)
colnames(perm_df) <- c("run1", "run2", "run3", "run4", "run5")

#Get the average q value from all of the permutations
perm_df <- perm_df %>% mutate(average_perm_q = rowMeans(.))

#Obtain q0 and q1 vectors
q0 <- perm_df$average_perm_q
q1 <- results_sorted$q

null.q <- q0
obs.q <- q1

#Cap q-values at min.q for drawing purposes
min.q <- 1e-16
obs.q[obs.q < min.q]  <- min.q
null.q[null.q < min.q] <- min.q

#Plot qq plot
n.test <- 57434
qqplot(-log10(null.q), -log10(obs.q), xlab = "Permuted Log10(q-value)", ylab = "Log10(q-value) caQTL")
abline(a=0, b=1)


#This function returns the P-value threshold corresponding to FDR=alpha.
getFDR <-
  function(q1, q0, alpha=0.05, z=NULL, subset=NULL){
    if(is.null(z)){
      a=0
      for(itr in 1:10){
        a=getFDR(q1,q0,alpha,rev(a+0:100/100^itr),subset)
      }
      a
    }else{
      if(!is.null(subset)){
        q1=q1[subset]
        q0=q0[subset]
      }
      q1=q1[!is.na(q1)]
      q0=q0[!is.na(q0)]
      x=NULL;
      for(i in z){
        x=c(x,sum(q0<i)/length(q0)/(sum(q1<i)/length(q1)))
      };
      max(c(0,z[x<alpha]),na.rm=T)
    }
  }

#Mark significant caQTLs (5% FDR)
#True = significant QTLs
flag_fdr05 = q1 < getFDR(q1, q0, 0.05) 
table(flag_fdr05)

#Attach flag_fdr05 column
results_sorted <- results_sorted %>% mutate(keep = flag_fdr05)

#Filter for rows passing FDR cutoff
FDR_results_endothelial <- results_sorted %>% filter(keep == TRUE)

#Write results
write.table(FDR_results_endothelial, file = "RASQUAL_results_endothelial_FDR_0.05.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
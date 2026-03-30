library(admixtools)
library(dplyr)
library(tidyr)
library(stringr)
library(stringr)
library(patchwork)


setwd("~/Library/CloudStorage/Box-Box/projects/dingo/new_dingo/data/")

getMergedDf <- function(suffix) {
  popdrop <- read.csv(paste0("popdrop_", suffix, ".csv"))
  weights <- read.csv(paste0("weights_", suffix, ".csv"))
  
  # Make weights wider
  weights_wider <- weights %>%
    pivot_wider(
      id_cols = c(target, source, right),
      names_from = left,
      values_from = c(weight, se, z),
      names_sep = "_") 
  
  # Ignore untested patterns
  popdrop_pat0 <- popdrop %>% subset(pat == 0)
  
  # Merge the two dataframes
  merged_df <- merge(popdrop_pat0, weights_wider, by=c("target", "source", "right")) %>%
    mutate(num_sources = str_count(source, "\\|") + 1,
           num_right = str_count(right, "\\|") + 1) %>%
    subset(feasible==TRUE)
  
  return(merged_df)
}


rankModels <- function(df) {
  all_z_scores_above_threshold <- function(row, z_columns) {
    #print(row)
    #print(z_columns)
    sources <- str_split(row["source"], "\\|")[[1]]             # Extract sources
    z_cols <- paste0("z_", sources)                            # Create corresponding z_ column names
    if (all(z_cols %in% z_columns)) {                          # Ensure all z_cols exist in the dataframe
      all(as.numeric(row[z_cols]) > 2)                         # Check if all relevant z-scores are > 3
    } else {
      FALSE                                                    # Return FALSE if any z_cols are missing
    }
  }
  
 sum_z_scores_above_threshold <- function(row, z_columns) {
   #print(row)
   #print(z_columns)
   sources <- str_split(row["source"], "\\|")[[1]]             # Extract sources
   z_cols <- paste0("z_", sources)                            # Create corresponding z_ column names
   if (all(z_cols %in% z_columns)) {                          # Ensure all z_cols exist in the dataframe
     sum(as.numeric(row[z_cols]) > 2)                         # Check if all relevant z-scores are > 3
   } else {
     NA                                                    # Return FALSE if any z_cols are missing
   }
 }
  
  ranked_df <- df %>%
    # Group by target
    group_by(target) %>%
    rowwise() %>%                                            # Process each row independently
    mutate(
      all_above_3 = all_z_scores_above_threshold(cur_data(), names(dplyr::select(cur_data(), starts_with("z_")))),
      sum_above_3 = sum_z_scores_above_threshold(cur_data(), names(dplyr::select(cur_data(), starts_with("z_"))))
    ) %>%
    #ungroup() %>%   
    # Arrange models within each target based on criteria
    group_by(target) %>%
    arrange(
      target,
      num_sources,                      # Fewer sources preferred
      desc(sum_above_3),              # Models with all z-scores > 3 preferred
      desc(num_right),
      desc(p)                                 # higher p-values are preferred as a tie-breaker
    ) %>%
    # Add rank column
    mutate(rank = row_number())
  
  
  return(ranked_df)
}



qpadm_hwd_pap <- getMergedDf("ngsd_chapter4_all_snps_v3") %>% 
  filter(!str_detect(source, "Highland")) %>%
  arrange(target, num_sources, desc(p)) %>%
  mutate(
    valid := case_when(
      p >= 0.01 ~ T,
      .default = F)) %>%
  dplyr::select(target, source, p, Dingo_Ancient_Nullarbor, HighlandWildDog, Papua_Indonesia_518bp, 
                se_Dingo_Ancient_Nullarbor, se_HighlandWildDog, se_Papua_Indonesia_518bp, z_Dingo_Ancient_Nullarbor,
                z_HighlandWildDog, z_Papua_Indonesia_518bp, num_sources, valid)

write.csv(qpadm_hwd_pap, "qpadm_hwd_pap_models.csv", row.names = F, quote = F)

the_rejects <- qpadm_hwd_pap %>%
  filter(!target %in% qpadm_hwd_pap_ranked_best$target)

qpadm_hwd_pap_ranked <- rankModels(getMergedDf("ngsd_chapter4_all_snps_v3") %>% filter(p > 0.01, !str_detect(source, "Highland")))
qpadm_hwd_pap_ranked_best <- qpadm_hwd_pap_ranked %>% filter(rank == 1)
write.csv(qpadm_hwd_pap_ranked_best, "qpadm_hwd_pap_ngsd_ranked_best.csv", quote = F, row.names = F)



qpadm_hwd_pap_ngsd <- getMergedDf("ngsd_chapter4") %>% 
  filter(p > 0.01) 
qpadm_hwd_pap_ngsd_ranked <- rankModels(qpadm_hwd_pap_ngsd)



qpadm_hwd_pap <- getMergedDf("ngsd_chapter4_v2") %>% 
  filter(p > 0.01)
qpadm_hwd_pap_ranked <- rankModels(qpadm_hwd_pap)

qpadm_vd <- getMergedDf("results_vd_more")

library(admixtools)
library(dplyr)
library(tidyr)

# # Capture command-line arguments
# args <- commandArgs(trailingOnly = TRUE)
# if (length(args) == 0) {
#   stop("No arguments provided. Please specify arguments.")
# }

# # Print all arguments
# cat("Arguments passed:", args, "\n")

#target <- args[1]
targets <- read.csv("villagedogs_targets.list", header=FALSE, col.names=c("ID"))$ID
#targets <- c("BorderCollie_BC0480")
#targets <- c("Dingo_Alpine_722g", "Dingo_Alpine_D01", "Dingo_Alpine_D05", "Dingo_FraserIsland_D02", 
#  "Dingo_Gibson_Desert_D08", "Dingo_Kimberley_D09", "Dingo_Northern_Queensland_D04", 
#  "Dingo_Northern_Queensland_D10", "Dingo_Northern_Territory_D06", "Dingo_Simpson_Desert_D03", "Dingo_Alpine_Cooinda",
#  "Dingo_Desert_Sandy", "Dingo_Alpine_D07")


setwd("~/Library/CloudStorage/Box-Box/projects/dingo/new_dingo/qpadm")

prefix <- "Sahul_canids_600K_Tvs"

#left_right <- c("CockerSpaniel", "Nullarbor", "Curracurrang", "NGSD")
left_right <- c("GermanShepherdDog", "Dingo_Ancient_Nullarbor", "Zhokhov9500BP.CGG6", "Germany_7k", "Ireland_Neolithic.Newgrange", "Russia_Baikal_7k", "Karelia_Veretye", "America_pool", "Sweden_PWC", "Iran_ChL")
right_fix <- c("CoyoteCalifornia") # Karelia_Veretye, 
#right_fix <- c("CoyoteCalifornia", "AndeanFox", "Baikal", "Basenji", "ChineseIndigenousDog", "IndigenousDogNigeria")



#all_combinations <- unlist(lapply(seq_along(left_right), function(x) combn(left_right, x, simplify = FALSE)), recursive = FALSE)
all_combinations <- unlist(
  lapply(seq_len(3), function(x) combn(left_right, x, simplify = FALSE)), 
  recursive = FALSE
)

expected_columns <- c(
  "pat", "wt", "dof", "chisq", "p", "f4rank",
  left_right,
  "feasible", "best", "dofdiff", "chisqdiff", "p_nested",
  "target", "right", "source"
)

# Create an empty dataframe with these columns, filled with NAs
empty_popdrop_df <- as.data.frame(matrix(NA, nrow = 0, ncol = length(expected_columns)))
colnames(empty_popdrop_df) <- expected_columns

# Function to ensure each popdrop result has all expected columns
ensure_columns <- function(df, columns) {
  missing_cols <- setdiff(columns, colnames(df))
  for (col in missing_cols) {
    df[[col]] <- NA
  }
  df <- df[, columns, drop = FALSE]
  return(df)
}

qpadm_function <- function(prefix, target, left_list, right_list, allsnps=TRUE) {
  print(paste("Running qpadm with left:", paste(left_list, collapse = ", "), "and right:", paste(right_list, collapse = ", ")))
  result <- qpadm(prefix, left = left_list, right = right_list, target = target, 
                  auto_only = FALSE,  allsnps = allsnps, constrained = FALSE, fudge_twice=TRUE)
  print(result)
  weights <- result$weights %>%
    mutate(target = target,
           right = paste(right_list, collapse = "|"),
           source = paste(left_list, collapse = "|"))
  popdrop <- result$popdrop %>%
    mutate(target = target, 
           right = paste(right_list, collapse = "|"),
           source = paste(left_list, collapse = "|"))  
  return(list(weights = weights, popdrop = popdrop))
}


qpadm_rotate_cooler <- function(prefix, target, left_right, right_fix, allsnps=TRUE) {
  weights_list <- list()
  popdrop_list <- list()
  for (comb in all_combinations) {
    # Items to be in the left list
    modified_left <- comb
    
    # Remaining items to be added to right_fix
    remaining_items <- setdiff(left_right, comb)
    
    # New right_fix list
    modified_right_fix <- c(right_fix, remaining_items)
    
    # Run the qpadm function with the modified lists
    curr_result <- qpadm_function(prefix, target, modified_left, modified_right_fix, allsnps=allsnps)
    weights_list <- append(weights_list, list(curr_result$weights))
    popdrop_list <- append(popdrop_list, list(curr_result$popdrop))


    # curr_result <- qpadm_function(prefix, target, modified_left, right_fix, allsnps=allsnps)
    # weights_list <- append(weights_list, list(curr_result$weights))
    # popdrop_list <- append(popdrop_list, list(curr_result$popdrop))
  }
  return(list(weights = weights_list, popdrop = popdrop_list))
}

#rotate_cooler_results <- qpadm_rotate_cooler(prefix, target, left_right, right_fix)
#targets <- c("Dingo_Alpine_D01", "Dingo_Alpine_D05", "Dingo_Alpine_D07", "Dingo_Alpine_722g",
#  "Dingo_FraserIsland_D02", "Dingo_FraserIsland_E56", "Dingo_FraserIsland_Y07", "Dingo_FraserIsland_Y47",
#  "Dingo_Gibson_Desert_D08", "Dingo_Kimberley_D09", "Dingo_Northern_Queensland_D04", "Dingo_Northern_Queensland_D10",
#  "Dingo_Northern_Territory_D06", "Dingo_Simpson_Desert_D03")

# Initialize lists to collect all weights and popdrop results
all_weights <- list()
all_popdrop <- list()

# Loop over each target sample
for (target in targets) {
  print(paste("Processing target:", target))
  
  # Run qpadm_rotate_cooler for the current target
  rotate_cooler_results <- qpadm_rotate_cooler(prefix, target, left_right, right_fix, allsnps=TRUE)
  popdrop_results <- lapply(rotate_cooler_results$popdrop, ensure_columns, expected_columns)
  # Append the results to the respective lists
  all_weights <- append(all_weights, rotate_cooler_results$weights)
  all_popdrop <- append(all_popdrop, popdrop_results)
}

# Convert the lists to data frames
all_weights_df <- do.call(rbind, all_weights)
all_popdrop_df <- do.call(rbind, all_popdrop)

# Save the results to files
write.csv(all_weights_df, paste0("weights_results_vd_more.csv"), row.names = FALSE)
write.csv(all_popdrop_df, paste0("popdrop_results_vd_more.csv"), row.names = FALSE)
print("Finito WGS!")

# prefix <- "souilmi2024"
# # Initialize lists to collect all weights and popdrop results
# all_weights <- list()
# all_popdrop <- list()

# # Loop over each target sample
# #for (target in targets) {  
#   # Run qpadm_rotate_cooler for the current target
#   rotate_cooler_results <- qpadm_rotate_cooler(prefix, target, left_right, right_fix, allsnps=FALSE)
#   popdrop_results <- lapply(rotate_cooler_results$popdrop, ensure_columns, expected_columns)
#   # Append the results to the respective lists
#   all_weights <- append(all_weights, rotate_cooler_results$weights)
#   all_popdrop <- append(all_popdrop, popdrop_results)
# #}

# # Convert the lists to data frames
# all_weights_df <- do.call(rbind, all_weights)
# all_popdrop_df <- do.call(rbind, all_popdrop)

# # Save the results to files
# write.csv(all_weights_df, paste0("weights_results_wgs_", target, ".csv"), row.names = FALSE)
# write.csv(all_popdrop_df, paste0("popdrop_results_wgs_", target, ".csv"), row.names = FALSE)
# print("Finito!")



# for (target in targets) {
#   left <- c("GermanShepherdDog", "Dingo_Ancient_Curracurrang")
#   right <- c("CoyoteCalifornia", "Zhokhov9500BP.CGG6", "Germany_7k", "Ireland_Neolithic.Newgrange", "Russia_Baikal_7k")
#   result <- qpadm(prefix, left = left, right = right, target = target, auto_only = FALSE,  allsnps = FALSE, constrained = FALSE, fudge_twice=TRUE)
#   print(result)
#   left <- c("Dingo_Ancient_Curracurrang")
#   right <- c("CoyoteCalifornia", "Zhokhov9500BP.CGG6", "Germany_7k", "Ireland_Neolithic.Newgrange", "Russia_Baikal_7k", "GermanShepherdDog")
#   result <- qpadm(prefix, left = left, right = right, target = target, auto_only = FALSE,  allsnps = FALSE, constrained = FALSE, fudge_twice=TRUE)
#   print(result)
# }

library(admixtools)
library(dplyr)
library(tidyr)
library(stringr)
library(stringr)
library(patchwork)


setwd("~/Library/CloudStorage/Box-Box/projects/dingo/heritage_paper/qpadm_revisions/")

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



write.csv(rankModels(getMergedDf("snp_lab") %>% filter(p > 0.01)) %>% filter(rank==1), "~/Library/CloudStorage/Box-Box/projects/dingo/heritage_paper/publication_plots/data/qpadm_snp_lab.csv",
          quote = F, row.names = F)
write.csv(rankModels(getMergedDf("snp_bordercollie") %>% filter(p > 0.01)) %>% filter(rank==1), "~/Library/CloudStorage/Box-Box/projects/dingo/heritage_paper/publication_plots/data/qpadm_snp_bordercollie.csv",
          quote = F, row.names = F)
write.csv(rankModels(getMergedDf("snp_scottterrier") %>% filter(p > 0.01)) %>% filter(rank==1), "~/Library/CloudStorage/Box-Box/projects/dingo/heritage_paper/publication_plots/data/qpadm_snp_scottterrier.csv",
          quote = F, row.names = F)

write.csv(rankModels(getMergedDf("wgs_lab") %>% filter(p > 0.01)) %>% filter(rank==1), "~/Library/CloudStorage/Box-Box/projects/dingo/heritage_paper/publication_plots/data/qpadm_wgs_lab.csv",
          quote = F, row.names = F)
write.csv(rankModels(getMergedDf("wgs_collie") %>% filter(p > 0.01)) %>% filter(rank==1), "~/Library/CloudStorage/Box-Box/projects/dingo/heritage_paper/publication_plots/data/qpadm_wgs_collie.csv",
          quote = F, row.names = F)
write.csv(rankModels(getMergedDf("wgs_scottterrier") %>% filter(p > 0.01)) %>% filter(rank==1), "~/Library/CloudStorage/Box-Box/projects/dingo/heritage_paper/publication_plots/data/qpadm_wgs_terrier.csv",
          quote = F, row.names = F)

df_sim <- getMergedDf("simulated_wgs_samples") %>%
  filter(num_sources > 1)

df_sim_final <- rankModels(df_sim %>% filter(p > 0.01)) %>% filter(rank == 1)
write.csv(df_sim_final, "~/Library/CloudStorage/Box-Box/projects/dingo/heritage_paper/publication_plots/data/qpadm_simulated_wo_single_source.csv",
          quote = F, row.names = F)


setwd("~/Library/CloudStorage/Box-Box/projects/dingo/scarsbrook_et_al_2025/qpadm")
df10 <- getMergedDf("results_scarsbrook_ngsd")
df11 <- getMergedDf("results_shyam_ngsd")
df12 <- getMergedDf("results_scarsbrook_kangaroodog")
df13 <- getMergedDf("results_ancient_scarsbrook")
length(unique(df13$target))


df <- getMergedDf("results_targets1") %>%
  filter(p > 0.01)
df2 <- getMergedDf("results_ngsd") %>%
  filter(p > 0.01)
df3 <- getMergedDf("results_ngsd_modern") %>%
  filter(p > 0.01)
df4 <- getMergedDf("results_ngsd_modern_gatton") %>%
  filter(p > 0.01)
df5 <- getMergedDf("results_snp_array_snps_g_n_c") %>%
  filter(p > 0.01)
df6 <- getMergedDf("results_vd_g_n_p") %>%
  filter(p > 0.01)
df7 <- getMergedDf("results_nhwd_g_n_p") %>%
  filter(p > 0.01)
df8 <- getMergedDf("results_vd_more") %>%
  filter(p > 0.01)
df9 <- getMergedDf("results_wgs_village_dogs")

########
setwd("~/Library/CloudStorage/Box-Box/projects/dingo/heritage_paper/qpadm_revisions")
#### Zhang
snp_w_ec_zhang <- rankModels(getMergedDf("snp_zhang_w_ec") %>% filter(p > 0.01)) %>% filter(rank == 1)
write.csv(snp_w_ec_zhang, 
          "~/Library/CloudStorage/Box-Box/projects/dingo/heritage_paper/publication_plots/data/qpadm_snp_array_zhang_w_ec.csv",
          quote=F,
          row.names = F)
snp_only_nulla_zhang <- rankModels(getMergedDf("snp_zhang_only_nulla") %>% filter(p > 0.01)) %>% filter(rank == 1)
write.csv(snp_only_nulla_zhang, 
          "~/Library/CloudStorage/Box-Box/projects/dingo/heritage_paper/publication_plots/data/qpadm_snp_array_zhang_only_nulla.csv",
          quote=F,
          row.names = F)

snp_only_curra_zhang <- rankModels(getMergedDf("snp_zhang_only_curra") %>% filter(p > 0.01)) %>% filter(rank == 1)
write.csv(snp_only_curra_zhang, 
          "~/Library/CloudStorage/Box-Box/projects/dingo/heritage_paper/publication_plots/data/qpadm_snp_array_zhang_only_curra.csv",
          quote=F,
          row.names = F)



snp_only_nulla <- getMergedDf("only_nulla")
snp_only_curra <- getMergedDf("only_curra")
wgs_only_nulla <- getMergedDf("wgs_only_nulla")
wgs_only_curra <- getMergedDf("wgs_only_curra")
snp_w_ec <- getMergedDf("w_ec")
wgs_w_ec <- getMergedDf("wgs_w_ec")


snp_w_ec_filt <- snp_w_ec %>% filter(p > 0.01)
wgs_w_ec_filt <- wgs_w_ec %>% filter(p > 0.01)
snp_only_nulla_filt <- snp_only_nulla %>% filter(p > 0.01)
wgs_only_nulla_filt <- wgs_only_nulla %>% filter(p > 0.01)
snp_only_curra_filt <- snp_only_curra %>% filter(p > 0.01)
wgs_only_curra_filt <- wgs_only_curra %>% filter(p > 0.01)

snp_w_ec_filt_ranked <- rankModels(snp_w_ec_filt)
wgs_w_ec_filt_ranked <- rankModels(wgs_w_ec_filt)
snp_only_nulla_filt_ranked <- rankModels(snp_only_nulla_filt)
wgs_only_nulla_filt_ranked <- rankModels(wgs_only_nulla_filt)
snp_only_curra_filt_ranked <- rankModels(snp_only_curra_filt)
wgs_only_curra_filt_ranked <- rankModels(wgs_only_curra_filt)

snp_w_ec_filt_ranked_top <- snp_w_ec_filt_ranked %>% filter(rank == 1)
wgs_w_ec_filt_ranked_top <- wgs_w_ec_filt_ranked %>% filter(rank == 1)
snp_only_nulla_filt_ranked_top <- snp_only_nulla_filt_ranked %>% filter(rank == 1)
wgs_only_nulla_filt_ranked_top <- wgs_only_nulla_filt_ranked %>% filter(rank == 1)
snp_only_curra_filt_ranked_top <- snp_only_curra_filt_ranked %>% filter(rank == 1)
wgs_only_curra_filt_ranked_top <- wgs_only_curra_filt_ranked %>% filter(rank == 1)

###### Scarsbrook
snp_scarsbrook <- getMergedDf("snp_modern_scarsbrook")
wgs_scarsbrook <- getMergedDf("wgs_modern_scarsbrook")
snp_only_nulla_scarsbrook <- getMergedDf("snp_modern_scarsbrook_only_nulla")
snp_only_curra_scarsbrook <- getMergedDf("snp_modern_scarsbrook_only_curra")
wgs_only_nulla_scarsbrook <- getMergedDf("wgs_modern_scarsbrook_only_nulla")
wgs_only_curra_scarsbrook <- getMergedDf("wgs_modern_scarsbrook_only_curra")
snp_w_ec_scarsbrook <- getMergedDf("snp_modern_scarsbrook_w_ec")
wgs_w_ec_scarsbrook <- getMergedDf("wgs_modern_scarsbrook_w_ec")

snp_scarsbrook_filt <- snp_scarsbrook %>% filter(p > 0.01)
wgs_scarsbrook_filt <- wgs_scarsbrook %>% filter(p > 0.01)
snp_w_ec_filt_scarsbrook <- snp_w_ec_scarsbrook %>% filter(p > 0.01)
wgs_w_ec_filt_scarsbrook <- wgs_w_ec_scarsbrook %>% filter(p > 0.01)
snp_only_nulla_filt_scarsbrook <- snp_only_nulla_scarsbrook %>% filter(p > 0.01)
wgs_only_nulla_filt_scarsbrook <- wgs_only_nulla_scarsbrook %>% filter(p > 0.01)
snp_only_curra_filt_scarsbrook <- snp_only_curra_scarsbrook %>% filter(p > 0.01)
wgs_only_curra_filt_scarsbrook <- wgs_only_curra_scarsbrook %>% filter(p > 0.01)

snp_scarsbrook_filt_ranked <- rankModels(snp_scarsbrook_filt)
wgs_scarsbrook_filt_ranked <- rankModels(wgs_scarsbrook_filt)
snp_w_ec_filt_ranked_scarsbrook <- rankModels(snp_w_ec_filt_scarsbrook)
wgs_w_ec_filt_ranked_scarsbrook <- rankModels(wgs_w_ec_filt_scarsbrook)
snp_only_nulla_filt_ranked_scarsbrook <- rankModels(snp_only_nulla_filt_scarsbrook)
wgs_only_nulla_filt_ranked_scarsbrook <- rankModels(wgs_only_nulla_filt_scarsbrook)
snp_only_curra_filt_ranked_scarsbrook <- rankModels(snp_only_curra_filt_scarsbrook)
wgs_only_curra_filt_ranked_scarsbrook <- rankModels(wgs_only_curra_filt_scarsbrook)

snp_scarsbrook_filt_ranked_top <- snp_scarsbrook_filt_ranked %>% filter(rank == 1)
wgs_scarsbrook_filt_ranked_top <- wgs_scarsbrook_filt_ranked %>% filter(rank == 1)
snp_w_ec_filt_ranked_top_scarsbrook <- snp_w_ec_filt_ranked_scarsbrook %>% filter(rank == 1)
wgs_w_ec_filt_ranked_top_scarsbrook <- wgs_w_ec_filt_ranked_scarsbrook %>% filter(rank == 1)
snp_only_nulla_filt_ranked_top_scarsbrook <- snp_only_nulla_filt_ranked_scarsbrook %>% filter(rank == 1)
wgs_only_nulla_filt_ranked_top_scarsbrook <- wgs_only_nulla_filt_ranked_scarsbrook %>% filter(rank == 1)
snp_only_curra_filt_ranked_top_scarsbrook <- snp_only_curra_filt_ranked_scarsbrook %>% filter(rank == 1)
wgs_only_curra_filt_ranked_top_scarsbrook <- wgs_only_curra_filt_ranked_scarsbrook %>% filter(rank == 1)


##### Merge Scarsbrook and Rest...
snp_w_ec_filt_ranked_top_merged <- rbind(snp_w_ec_filt_ranked_top, snp_w_ec_filt_ranked_top_scarsbrook)
wgs_w_ec_filt_ranked_top_merged <- rbind(wgs_w_ec_filt_ranked_top, wgs_w_ec_filt_ranked_top_scarsbrook)
snp_only_nulla_filt_ranked_top_merged <- rbind(snp_only_nulla_filt_ranked_top, snp_only_nulla_filt_ranked_top_scarsbrook)
wgs_only_nulla_filt_ranked_top_merged <- rbind(wgs_only_nulla_filt_ranked_top, wgs_only_nulla_filt_ranked_top_scarsbrook)
snp_only_curra_filt_ranked_top_merged <- rbind(snp_only_curra_filt_ranked_top, snp_only_curra_filt_ranked_top_scarsbrook)
wgs_only_curra_filt_ranked_top_merged <- rbind(wgs_only_curra_filt_ranked_top, wgs_only_curra_filt_ranked_top_scarsbrook)


write.csv(wgs_scarsbrook_filt_ranked %>% dplyr::select(target, source, right, p, GermanShepherdDog, Dingo_Ancient_Curracurrang, Dingo_Ancient_Nullarbor, se_GermanShepherdDog, se_Dingo_Ancient_Curracurrang, se_Dingo_Ancient_Nullarbor, z_GermanShepherdDog, z_Dingo_Ancient_Curracurrang, z_Dingo_Ancient_Nullarbor ,num_sources, rank),
    "supp_table_sacrsbrook.csv", 
    quote=F,
    row.names = F
          )

write.csv(snp_scarsbrook_filt_ranked_top, 
          "~/Library/CloudStorage/Box-Box/projects/dingo/heritage_paper/publication_plots/data/qpadm_snp_array_scarsbrook.csv",
          quote=F,
          row.names = F)
write.csv(wgs_scarsbrook_filt_ranked_top, 
          "~/Library/CloudStorage/Box-Box/projects/dingo/heritage_paper/publication_plots/data/qpadm_wgs_array_scarsbrook.csv",
          quote=F,
          row.names = F)

write.csv(snp_w_ec_filt_ranked_top_merged, 
          "~/Library/CloudStorage/Box-Box/projects/dingo/heritage_paper/publication_plots/data/qpadm_snp_array_cockerspaniel.csv",
          quote=F,
          row.names = F)
write.csv(wgs_w_ec_filt_ranked_top_merged, 
          "~/Library/CloudStorage/Box-Box/projects/dingo/heritage_paper/publication_plots/data/qpadm_wgs_cockerspaniel.csv",
          quote=F,
          row.names = F)
write.csv(snp_only_nulla_filt_ranked_top_merged, 
          "~/Library/CloudStorage/Box-Box/projects/dingo/heritage_paper/publication_plots/data/qpadm_snp_array_only_nulla.csv",
          quote=F,
          row.names = F)
write.csv(wgs_only_nulla_filt_ranked_top_merged, 
          "~/Library/CloudStorage/Box-Box/projects/dingo/heritage_paper/publication_plots/data/qpadm_wgs_only_nulla.csv",
          quote=F,
          row.names = F)
write.csv(snp_only_curra_filt_ranked_top_merged, 
          "~/Library/CloudStorage/Box-Box/projects/dingo/heritage_paper/publication_plots/data/qpadm_snp_array_only_curra.csv",
          quote=F,
          row.names = F)
write.csv(wgs_only_curra_filt_ranked_top_merged, 
          "~/Library/CloudStorage/Box-Box/projects/dingo/heritage_paper/publication_plots/data/qpadm_wgs_only_curra.csv",
          quote=F,
          row.names = F)

df13_p_filt <- df13 %>% filter(p > 0.01)
length(unique(df13_p_filt$target))

df13_p_filt_ranked <- rankModels(df13_p_filt)
write.csv(df13_p_filt_ranked, "results_ancient_scarsbrook_ranked.csv", quote = F, row.names = F)
write.csv(df13_p_filt_ranked %>% filter(rank == 1), "results_ancient_scarsbrook_top.csv", quote = F, row.names = F)

df14 <- getMergedDf("results_ancient_scarsbrook_highcov_simple")
length(unique(df14$target))

df15 <- getMergedDf("results_ancient_scarsbrook_highcov_simple_ancdog")
length(unique(df15$target))

df14_p_filt <- df14 %>% filter(p > 0.01)
length(unique(df14_p_filt$target))

df15_p_filt <- df15 %>% filter(p > 0.01)
length(unique(df15_p_filt$target))

df16 <- getMergedDf("results_ancient_scarsbrook_highcov_ancdog_wpapua")
length(unique(df16$target))

df16_p_filt <- df16 %>% filter(p > 0.01)
length(unique(df16_p_filt$target))

df17 <- getMergedDf("results_ancient_scarsbrook_highcov_troubleshoot")
length(unique(df17$target))

df17_p_filt <- df17 %>% filter(p > 0.01)
length(unique(df17_p_filt$target))

df17_p_filt_ranked <- rankModels(df17_p_filt)

df3_ranked <- rankModels(df3)
df3_ranked %>% filter(rank == 1)

write.csv(df3_ranked, "ngsd_modern.csv", quote = F, row.names = F)

missing_samples <- df %>% 
  filter(target %in% c("Dingo_Ancient_Gatton_14", "Dingo_Ancient_VIC_8", "Dingo_Ancient_VIC_5", "Papua_Indonesia_518bp"))

write.csv(df, "targets1_qpadm.csv", quote = F, row.names = F )
write.csv(df2, "ngsd_qpadm.csv", quote=F, row.names=F)
write.csv(df3, "ngsd_ancient_modelling_modern.csv", quote=F, row.names=F)
write.csv(df4, "ngsd_ancient_modelling_modern_w_gatton_source.csv", quote=F, row.names=F)
write.csv(df8, "village_dogs_all_sources.csv", quote=F, row.names=F )

df_long <- df7 %>% 
  select(target, source, starts_with(c("weight", "se")), p) %>%
  pivot_longer(cols = starts_with(c("weight", "se")),
               names_to = c(".value", "source_label"),
               names_pattern = "(weight|se)_(.*)") %>%
  subset(!is.na(weight)) %>%
  group_by(target, source) %>% 
  arrange(desc(source_label), by.group = TRUE) %>%
  mutate(weight2 = cumsum(weight))


plot_one <- function(sample) {
  df_long %>% 
    subset(target==sample) %>%
    ggplot(aes(x = interaction(target, source), y = weight,fill = source_label)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_errorbar(aes(ymin = weight2 - se, ymax = weight2 + se), width=0.3, position = position_dodge(width = 0.5)) +
    geom_point(aes(y=weight2), position = position_dodge(width = 0.5)) +
    geom_text(aes(x=interaction(target, source), y=1.05, label=signif(p, digits = 2))) + 
    facet_grid(~target) + 
    coord_cartesian(ylim=c(-0.1, 1.1)) +
    scale_y_continuous(breaks=c(0,0.5,1)) + 
    guides(fill="none") +
    scale_fill_manual(values=manual_colours) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(), 
          panel.grid.major.x = element_blank()) 
}

# Make a list of plots for all targets
plot_list <- lapply(unique(df_long$target), plot_one)

legend_plot <- ggplot(data.frame(source_label = names(manual_colours), y = 1),
                      aes(x = y, y = y, fill = source_label)) +
  geom_col() +
  scale_fill_manual(values = manual_colours, name = "Source")

# Extract legend as grob
legend_only <- cowplot::get_legend(legend_plot)

# Wrap it so patchwork can use it
legend_panel <- patchwork::wrap_elements(legend_only)

# Combine your plots + legend
final_plot <- wrap_plots(
  plotlist = c(plot_list), 
  ncol = 10
)
final_plot

ggsave("surbakti_snps_dingo_models.png", final_plot, dpi=500, height=10, width=20)

manual_colours <- c("GermanShepherdDog"="#ef476f", 
                    "NewGuineaSingingDog"="#f78c6b",
                    "Papua_Indonesia_518bp"= "#ffd166",
                    "Dingo_Ancient_Gatton_13"="#06d6a0",
                    "Dingo_Ancient_Curracurrang"="#118ab2",
                    "Dingo_Ancient_Nullarbor"="#073b4c")



legend_df <- data.frame(source_label = names(manual_colours))

legend_plot <- ggplot(legend_df, aes(x = 1, y = source_label, color = source_label)) +
  geom_point(size = 5) +
  scale_color_manual(values = manual_colours, name = "Source") +
  theme_void() +
  theme(legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 6)))




df6
head(df6) 
shannon_summary <- df8 %>%
  mutate(pop = str_remove(str_replace(target, "^VD", ""), "_[^_]+$")) %>%
  mutate(across(everything(), ~replace_na(., 0))) %>%
  group_by(target, pop) %>%
  mutate(ancient_canid_ancestry=weight_Dingo_Ancient_Nullarbor+weight_Papua_Indonesia_518bp) %>%
  summarise(dog=mean(weight_GermanShepherdDog),
            canid=mean(ancient_canid_ancestry)) %>%
  group_by(pop) %>%
  summarise(pop_dog = mean(dog),
            pop_canid = mean(canid),
            num_samples = n(),
            se_dog = sd(dog),
            se_canid = sd(canid)) %>%
  mutate(sum = pop_dog+pop_canid)

df8_ranked <- rankModels(df8)
write.csv(df8_ranked, "ranked_df.csv", quote = F, row.names = F)


df7_ranked <- rankModels(df7)
write.csv(df7_ranked, "ranked_df.csv", quote = F, row.names = F)


best_model_per_sample <- df8_ranked %>% filter(rank == 1)
best_model_per_sample_sur <- df7_ranked %>% filter(rank == 1)

write.csv(best_model_per_sample, "../shannon_boyoko_et_al/best_models.csv", quote = F, row.names = F)
write.csv(best_model_per_sample_sur, "best_models.csv", quote = F, row.names = F)



write.csv(shannon_summary, "shannon_summary.csv", quote = F, row.names = F)


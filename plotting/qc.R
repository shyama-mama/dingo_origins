pacman::p_load(
  dplyr,
  ggplot2,
  stringr,
  rnaturalearth,
  rnaturalearthdata,
  tidyr
)


setwd("~/Library/CloudStorage/Box-Box/projects/dingo/new_dingo/data")

masked_missingness <- fread("merge_masked_all_v2_missing.imiss") %>%
  janitor::clean_names()

mosaic_dingo <- fread("~/Library/CloudStorage/Box-Box/projects/dingo/MOSAIC/results_batch260114/mosaic_wgs_dingoes.csv") %>%
  janitor::clean_names()
mosaic_sahul <- fread("~/Library/CloudStorage/Box-Box/projects/dingo/MOSAIC/results_sahul/mosaic_wgs.csv") %>%
  janitor::clean_names()

mosaic_all_wgs <- rbind(mosaic_dingo, mosaic_sahul)

metadata <- fread("~/Library/CloudStorage/Box-Box/projects/dingo/new_dingo/metadata/metadata_for_wgs_samples.tsv") %>%
  janitor::clean_names()

mosaic_missingness <- merge(mosaic_all_wgs, masked_missingness, by.x="sample", by.y="fid") %>%
  left_join(., metadata, by="sample")


##### Check the missingness of samples
ggplotly(mosaic_missingness %>%
  ggplot() +
  geom_abline(intercept = 0) +
  geom_hline(yintercept = 0.50) + 
  geom_point(aes(x=dingo_ancestry, y=1-f_miss, label=sample, colour=population_label)) + 
  theme_bw() +
  theme(
    panel.grid = element_blank()
  ))


###### Relatedness 
king_relatedness <- fread("merge_unmasked.kin0") %>%
  janitor::clean_names()

related_individuals <- king_relatedness %>%
  filter(kinship > 0.0884)

samples_to_remove <- c("Dingo_Alpine_722g", 
                       "Dingo_Alpine_D07", 
                       "Dingo15",
                       "di1426",
                       "DO-CO-167564",
                       "DO-CO-167767",
                       "DO-CO-168535",
                       "ND1039",
                       "Dingo18",
                       "Dingo13",
                       "Dingo14",
                       "Dingo16",
                       "Dingo17",
                       "Dingo20",
                       "NewGuineaSingingDog",
                       "NewGuineaSingingDog07",
                       "NewGuineaSingingDog02",
                       "NewGuineaSingingDog03",
                       "NewGuineaSingingDog_BNGS",
                       "NewGuineaSingingDog16",
                       "NewGuineaSingingDog08",
                       "NewGuineaSingingDog14",
                       "NewGuineaSingingDog13",
                       "HighlandWildDog01")



related_individuals %>%
  filter(!fid1 %in% samples_to_remove, !fid2 %in% samples_to_remove)

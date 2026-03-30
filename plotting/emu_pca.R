pacman::p_load(
  dplyr,
  data.table,
  plotly,
  umap,
  tsne,
  scatterpie,
  rnaturalearth,
  rnaturalearthdata,
  sf,
  ggpubr
)


setwd("~/Library/CloudStorage/Box-Box/projects/dingo/new_dingo/data")


sample_order <- fread("merge_masked_all_v3_dingoes_unrelated_filt_geno0.8_maf0.01_thin1000.fam", header = F,
                      col.names = c("pop_id","sample_id","mum","dad","sex","case_control")) %>%
  dplyr::select(-c(mum, dad, sex, case_control)) 

eigenvals <- fread("merge_masked_all_v3_dingoes_unrelated_filt_geno0.8_maf0.01_thin1000_emu.eigenvals", header=F) %>%
  mutate(pct = (V1/sum(V1))*100)
eigenvecs <- fread("merge_masked_all_v3_dingoes_unrelated_filt_geno0.8_maf0.01_thin1000_emu.eigenvecs", header=F)

metadata <- fread("../metadata/metadata_for_wgs_samples.tsv") %>%
  janitor::clean_names() 

sample_order_metadata <- left_join(sample_order, metadata, by=join_by("pop_id"=="sample")) 

pca_results <- cbind(sample_order_metadata, eigenvecs)

pca_results 
eigenvals %>%
  mutate(pct = (V1/sum(V1))*100)

ggplotly(pca_results %>%
           ggplot() +
           geom_point(aes(x=V1, y=V2, label=pop_id, colour=population_label)))

ggplotly(pca_results %>%
           ggplot() +
           geom_point(aes(x=V2, y=V3, label=pop_id, colour=population_label)))

ggplotly(pca_results %>%
           ggplot() +
           geom_point(aes(x=V3, y=V4, label=pop_id, colour=population_label)))


pc1_pc2 <- ggplot() +
  geom_hline(yintercept = 0, alpha=0.2) +
  geom_vline(xintercept = 0, alpha=0.2) + 
  geom_point(data=pca_results %>% filter(age=="Modern", population_label != ""), aes(x=V1, y=V2, label=pop_id, colour=population_label), size=2) +
  geom_point(data=pca_results %>% filter(age=="Ancient"), aes(x=V1, y=V2, label=pop_id, fill=population_label), shape=24, size=2.5, alpha=0.7) +
  theme_bw() +
  guides(fill="none") + 
  xlab(paste0("PC1 (", round(eigenvals$pct[1],1), "%)"))+
  ylab(paste0("PC2 (", round(eigenvals$pct[2],1), "%)"))+
  scale_colour_manual("Population", values=pop_colours)+
  scale_fill_manual(values=pop_colours)+
  theme(
    panel.grid = element_blank()
  )

pc2_pc3 <- ggplot() +
  geom_hline(yintercept = 0, alpha=0.2) +
  geom_vline(xintercept = 0, alpha=0.2) + 
  geom_point(data=pca_results %>% filter(age=="Modern", population_label != ""), aes(x=V2, y=V3, label=pop_id, colour=population_label), size=2) +
  geom_point(data=pca_results %>% filter(age=="Ancient"), aes(x=V2, y=V3, label=pop_id, fill=population_label), shape=24, size=2.5, alpha=0.7) +
  theme_bw() +
  guides(fill="none") + 
  xlab(paste0("PC2 (", round(eigenvals$pct[2],1), "%)"))+
  ylab(paste0("PC3 (", round(eigenvals$pct[3],1), "%)"))+
  scale_colour_manual("Population", values=pop_colours)+
  scale_fill_manual(values=pop_colours)+
  theme(
    panel.grid = element_blank()
  )


with_ngsd <- ggarrange(pc1_pc2, pc2_pc3, ncol = 2, common.legend = T, legend = "right", labels = c("a", "b"))

######## without ngsd
sample_order_nongsd <- fread("merge_masked_all_v3_dingoes_unrelated_filt_geno0.8_maf0.01_thin1000_nongsd.fam", header = F,
                      col.names = c("pop_id","sample_id","mum","dad","sex","case_control")) %>%
  dplyr::select(-c(mum, dad, sex, case_control)) 

eigenvals_nongsd <- fread("merge_masked_all_v3_dingoes_unrelated_filt_geno0.8_maf0.01_thin1000_nongsd_emu.eigenvals", header=F) %>%
  mutate(pct = (V1/sum(V1))*100)
eigenvecs_nongsd <- fread("merge_masked_all_v3_dingoes_unrelated_filt_geno0.8_maf0.01_thin1000_nongsd_emu.eigenvecs", header=F)

sample_order_metadata_nongsd <- left_join(sample_order_nongsd, metadata, by=join_by("pop_id"=="sample")) 

pca_results_nongsd <- cbind(sample_order_metadata_nongsd, eigenvecs_nongsd)


pc1_pc2_nongsd <- ggplot() +
  geom_hline(yintercept = 0, alpha=0.2) +
  geom_vline(xintercept = 0, alpha=0.2) + 
  geom_point(data=pca_results_nongsd %>% filter(age=="Modern", population_label != ""), aes(x=V1, y=V2, label=pop_id, colour=population_label), size=2) +
  geom_point(data=pca_results_nongsd %>% filter(age=="Ancient"), aes(x=V1, y=V2, label=pop_id, fill=population_label), shape=24, size=2.5, alpha=0.7) +
  theme_bw() +
  guides(fill="none") + 
  xlab(paste0("PC1 (", round(eigenvals_nongsd$pct[1],1), "%)"))+
  ylab(paste0("PC2 (", round(eigenvals_nongsd$pct[2],1), "%)"))+
  scale_colour_manual("Population", values=pop_colours)+
  scale_fill_manual(values=pop_colours)+
  theme(
    panel.grid = element_blank()
  )

pc2_pc3_nongsd <- ggplot() +
  geom_hline(yintercept = 0, alpha=0.2) +
  geom_vline(xintercept = 0, alpha=0.2) + 
  geom_point(data=pca_results_nongsd %>% filter(age=="Modern", population_label != ""), aes(x=V2, y=V3, label=pop_id, colour=population_label), size=2) +
  geom_point(data=pca_results_nongsd %>% filter(age=="Ancient"), aes(x=V2, y=V3, label=pop_id, fill=population_label), shape=24, size=2.5, alpha=0.7) +
  theme_bw() +
  guides(fill="none") + 
  xlab(paste0("PC2 (", round(eigenvals_nongsd$pct[2],1), "%)"))+
  ylab(paste0("PC3 (", round(eigenvals_nongsd$pct[3],1), "%)"))+
  scale_colour_manual("Population", values=pop_colours)+
  scale_fill_manual(values=pop_colours)+
  theme(
    panel.grid = element_blank()
  )
no_ngsd <- ggarrange(pc1_pc2_nongsd, pc2_pc3_nongsd, ncol = 2, common.legend = T, legend = "right", labels = c("c", "d"))

combined <- ggarrange(with_ngsd, no_ngsd, nrow = 2, common.legend = T, legend = "right")
ggsave("../emu/final_pca_plots.png", combined, dpi=300, height=6, width=8)


pacman::p_load(
  dplyr,
  data.table,
  pheatmap,
  reshape,
  MASS,
  tidyr,
  ggplot2,
  stringr,
  grid
)

setwd("~/Library/CloudStorage/Box-Box/projects/dingo/new_dingo/data/")

outgroup_f3_all <- fread("f3_outgroup_v3.csv") %>% 
  dplyr::select(pop2, pop3,est,se,z,n) 

pca_samples <- fread("pca_samples.list", header=F, col.names=c("id"))$id

outgroup_f3_filter <- outgroup_f3_all %>%
  filter(n > 10000,
         pop2 %in% pca_samples,
         pop3 %in% pca_samples)



############
metadata <- fread("../metadata/metadata_for_wgs_samples.tsv") %>%
  janitor::clean_names()

preiod_sample <- metadata %>%
  filter(population_label != "",
         species %in% c("Dingo", "NGSD", "HWD"),
         sample %in% unique(c(outgroup_f3_filter$pop2, outgroup_f3_filter$pop3))) %>% 
  dplyr::select(sample, population_label, age)

outgroup_f3_filter_w_meta <- outgroup_f3_filter %>%
  filter(
    pop2 %in% preiod_sample$sample,
    pop3 %in% preiod_sample$sample)

cast_matrix <- cast(outgroup_f3_filter_w_meta %>% dplyr::select(pop2,pop3,est), pop2~pop3)
rownames(cast_matrix) <- cast_matrix$pop2
cast_matrix_as_matrix <- as.matrix(cast_matrix, rownames = "pop3")
diag(cast_matrix_as_matrix)=NA

new_pops <- data.frame("Population" = preiod_sample$population_label)
rownames(new_pops) <- preiod_sample$sample

age_sample <- data.frame("Age" = preiod_sample$age)
rownames(age_sample) <- preiod_sample$sample


pop_colours <- c("East"="#ffca3a",
               "North"="#8aaec4",
               "Alpine"="#8ac926",
               "Captive"="#ff924c",       
               "West"="#1982c4",
               "Mallee"="#ff595e",
               "K'gari"="#841c26",
               "Central"="#7F7F7F",
               "HWD"="grey",
               "NGSD"="grey30",
               "Ancient East"="#FFAE43",
               "Ancient West"="#0D4162",
               "Ancient South"="#68971D",
               "Ancient NG"="black")
colours <- list(
  Population = pop_colours,
  Age=c(
    "Ancient"="#b5c6e0",
    "Modern"="#ebf4f5"
  ))


f3_heatmap <- pheatmap(cast_matrix_as_matrix,
                       clustering_callback=function(...)hclust(dist(cast_matrix_as_matrix) ),
                       border_color=NA,
                       annotation_row = age_sample,
                       annotation_col= new_pops,
                       annotation_colors = colours,
                       annotation_names_col = TRUE)

f3_heatmap
# Define a color for the column names (e.g., "blue")
ggsave("../outgroup_f3/f3_heatmap_v3.png", f3_heatmap, dpi=300, height=13, width=18)


######### 
outgroup_f3_filter_w_meta <- outgroup_f3_filter_w_meta %>% mutate(inverse_Dstat=est^-1) %>%
  mutate(negDstat=1-est)

#### MDS plot ####
#extract cols for plot
MDS_f3 <- dplyr::select(outgroup_f3_filter_w_meta,pop2,pop3,negDstat)

#pivot wide
MDS <- pivot_wider(MDS_f3, names_from=pop2, values_from = negDstat, values_fill = NA)
#remove sample names column
MDS_matrix <- MDS[, -1]
#format as matrix
MDS_matrix <- dist(t(as.matrix(MDS_matrix)))

# caluculate MDS
mds1 <- cmdscale(MDS_matrix, k=2, eig=TRUE)

#set output as dataframe
mds_n <- as.data.frame(mds1[["points"]])
setDT(mds_n, keep.rownames = TRUE)
mds_n=mds_n[order(mds_n[,1]),]
mds_n_w_meta <- left_join(mds_n, preiod_sample, by=join_by("rn"=="sample"))

mds_plot <- ggplot(mds_n_w_meta, aes(x=V1, y=V2, label=rn)) +
  geom_point(aes(fill=population_label, shape=age), size=4, alpha=0.7) +
  scale_fill_manual("Population", values=pop_colours) +
  scale_shape_manual(values=c("Modern"=21, "Ancient"=24)) + 
  theme_bw() + 
  guides(fill=guide_legend(override.aes = list(shape = 21, color = "black", stroke = 1)))+
  geom_hline(yintercept = 0, color = "black", alpha=0.2, size=0.2) +
  geom_vline(xintercept = 0, color = "black", alpha=0.2, size=0.2) +
  labs(x="Dimension 1", y= "Dimension 2", 
       title = "Multidimensional Scaling Plot of 1-f3 Filtered") +
  theme(
    panel.grid = element_blank()
  )
mds_plot

ggsave("../outgroup_f3/mds_plot_v3.png", mds_plot, dpi=300, height=5, width=6)




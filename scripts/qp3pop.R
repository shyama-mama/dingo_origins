pacman::p_load(
  admixtools,
  ggplot2,
  dplyr,
  tidyr,
  tidyverse,
  pheatmap,
  reshape,
  stringr,
  MASS,
  data.table,
  rnaturalearth,
  plotly,
  ggrepel
)

setwd("~/Library/CloudStorage/Box-Box/projects/dingo/new_dingo/outgroup_f3")


data <- "merge_masked_all_f3"
modern_dingoes <- fread("modern_dingoes.list", header=F, col.names = "id")$id
sahul_dogs <- fread("sahul_dogs.list", header=F, col.names = "id")$id
ancient_dingoes_f3 <- fread("ancient_dingoes_for_f3.list", header=F, col.names = "id")$id
pops <- c(modern_dingoes, sahul_dogs, ancient_dingoes_f3)
### Outgroup f3
outgroup <- "CoyoteCalifornia"

f3_matrix <- qp3pop("merge_masked_all_f3", pop1=outgroup, pop2=pops, pop3=pops, 
                    verbose = TRUE, 
                    adjust_pseudohaploid = TRUE,
                    allsnps=TRUE,
                    auto_only=FALSE,
                    boot=FALSE,
                    outgroupmode=TRUE)


write.table(f3_matrix_vic, file="f3_outgroup_souilmi2024_bergstrom_plassias_shyam2025_transversions_for_qpadm_qpgraph_vic.csv", sep = ",", row.names=FALSE, quote=FALSE)

##### f3
f3_matrix <- fread("f3_outgroup_souilmi2024_bergstrom_plassias_shyam2025_transversions_for_qpadm_qpgraph.csv")

bad_samples <- c("NewGuineaSingingDog_BNGS1", "NewGuineaSingingDog_BNGS2", "Dingo_Ancient_Foreshore_697bp", "Dingo_Ancient_Foreshore_715bp", "Dingo_Ancient_Gatton_LP_14")

cast_matrix <- cast(f3_matrix %>% dplyr::select(pop2,pop3,est) %>% filter(!pop2 %in% bad_samples,
                                                                          !pop3 %in% bad_samples), pop2~pop3)
rownames(cast_matrix) <- cast_matrix$pop2
cast_matrix_as_matrix <- as.matrix(cast_matrix, rownames = "pop3")
diag(cast_matrix_as_matrix)=NA

f3_heatmap <- pheatmap(cast_matrix_as_matrix,
                       clustering_callback=function(...)hclust(dist(cast_matrix_as_matrix) ),
                       border_color=NA)
f3_heatmap
ggsave("outgroup_f3_heatmap.png", f3_heatmap, dpi=500, height=10, width=10)
f3_matrix <- f3_matrix %>% 
  filter(!pop2 %in% bad_samples,
         !pop3 %in% bad_samples) %>%
  mutate(negDstat=1-est)

#### MDS plot ####
#extract cols for plot
MDS_f3 <- dplyr::select(f3_matrix,pop2,pop3,negDstat)

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


ggplotly(ggplot(mds_n, aes(x=V1, y=V2, label=rn)) +
  geom_point(size=4, alpha=0.7))

###### mds with only dingoes and ngsd
dogs <- c("Dingo_Ancient_VIC_10", "Dingo_Ancient_VIC_05", "Dingo_Ancient_VIC_04", "GermanShepherdDog10", "GermanShepherdDog_DS051", "GermanShepherdDog_DS053", "GermanShepherdDog12")

f3_matrix_dingoes <- f3_matrix %>%
  filter(!pop2 %in% dogs,
         !pop3 %in% dogs)

MDS_f3_dingoes <- dplyr::select(f3_matrix_dingoes,pop2,pop3,negDstat)

#pivot wide
MDS <- pivot_wider(MDS_f3_dingoes, names_from=pop2, values_from = negDstat, values_fill = NA)
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

mds_n_metadata <- mds_n %>%
  mutate(type := case_when(
    str_detect(rn, "Ancient") ~ "Ancient",
    str_detect(rn, "Papua_Indonesia_518bp") ~ "Ancient",
    .default = "Modern"
  ),
  region := case_when(
    str_detect(rn, "Dingo_Northern_Territory_D06") ~ "Captive",
    str_detect(rn, "Cooinda") ~ "Captive",
    str_detect(rn, "Dingo_Alpine_D01") ~ "Captive",
    str_detect(rn, "Alpine") ~ "Alpine",
    str_detect(rn, "Curracurrang") ~ "Curracurrang",
    str_detect(rn, "Nullarbor") ~ "Nullarbor",
    str_detect(rn, "Gatton") ~ "Coastal NSW",
    str_detect(rn, "Ancient_NSW") ~ "Coastal NSW",
    str_detect(rn, "Dingo_Ancient_Gambier_724bp") ~ "Mt. Gambier",
    str_detect(rn, "Dingo_Ancient_VIC") ~ "Victoria", 
    str_detect(rn, "Dingo_FraserIsland") ~ "K'gari",
    str_detect(rn, "Dingo_Desert_Sandy") ~ "Desert",
    str_detect(rn, "Dingo_Gibson_Desert_D08") ~ "Desert",
    str_detect(rn, "Dingo_Kimberley_D09") ~ "Northern",
    str_detect(rn, "Dingo_Northern_Queensland") ~ "Northern",
    str_detect(rn, "Dingo_Simpson_Desert_D03") ~ "Desert",
    .default = "New Guinea Singing Dog"
  ))



# write.csv(mds_n_metadata, "mds_n_metadata.csv", quote=F, row.names = F)

mds_n_metadata_extra <- fread("mds_n_metadata_extra.csv")
ancient_modern_shape <- c("Ancient"=21, "Modern"=23)
dingo_ngsd_pops_colours <- c("New Guinea Singing Dog"="#000000", 
                             "Captive"="#E69F00", 
                             "K'gari"="#56B4E9", 
                             "Desert"="#009E73", 
                             "East"="#F0E442", 
                             "South"="#0072B2", 
                             "North"="#D55E00", 
                             "Nullarbor"="#CC79A7", 
                             "Alpine"="#DC3220")



mds_n_metadata_extra_pop <- mds_n_metadata_extra %>%
  mutate(population := case_when(
    region %in% c("Mt. Gambier", "Victoria") ~ "South",
    region %in% c("Coastal NSW", "Curracurrang") ~ "East",
    region %in% c("Northern") ~ "North",
    .default = region
  ),
  new = case_when(
    rn %in% c("Dingo_Ancient_Gatton_LP_13", "Dingo_Ancient_Gambier_724bp", "Dingo_Ancient_VIC_06", "Dingo_Ancient_VIC_08", "Papua_Indonesia_518bp", "Dingo_Ancient_NSW_LP_12") ~ "New",
    .default = "Published"
  ))

dingo_ngsd_pca <- ggplot(mds_n_metadata_extra_pop, aes(x=V1, y=V2, label=rn, fill=population, shape=type, colour=new, size=new)) +
  geom_point(stroke=1, alpha=0.8) +
  scale_shape_manual(name="Type", values=ancient_modern_shape) +
  scale_fill_manual(name="Population", values=dingo_ngsd_pops_colours) +
  scale_size_manual(name="Sample Origin", values = c("New"=5, "Published"=4)) + 
  scale_colour_manual(name="Sample Origin", values=c("New"="red", "Published"="black")) +
  #guides(fill=guide_legend(override.aes = list(shape = 21)), colour=guide_legend(override.aes = list(shape = 21))) +
  guides(fill="none", shape="none", colour="none", size="none") + 
  labs(x="MDS1", y="MDS2") + 
  theme_classic()
dingo_ngsd_pca
#ggplotly(dingo_ngsd_pca)

ggsave("dingo_ngsd_pca.png", dingo_ngsd_pca, dpi=500)


admixture_pc1 <- mds_n_metadata_extra_pop %>%
  filter(type == "Modern",
         region != "New Guinea Singing Dog") %>%
  ggplot(aes(x=V1, y=Admixture)) +
  geom_smooth(method="lm", linetype = "dashed", size=0.5, colour="black", alpha=0.3) +
  geom_point(aes(fill=population, shape=type), size=3, alpha=0.8) + 
  labs(x="MDS1", y="Dingo Ancestry") + 
  scale_shape_manual(name="Type", values=ancient_modern_shape) +
  scale_fill_manual(name="Population", values=dingo_ngsd_pops_colours) +
  scale_colour_manual(name="Sample Origin", values=c("New"="red", "Published"="black")) +
  guides(fill="none", shape="none", colour="none") + 
  #guides(fill=guide_legend(override.aes = list(shape = 21)), colour=guide_legend(override.aes = list(shape = 21))) +
  theme_classic()

admixture_pc2 <- mds_n_metadata_extra_pop %>%
  filter(type == "Modern",
         region != "New Guinea Singing Dog") %>%
  ggplot(aes(x=V2, y=Admixture)) +
  geom_smooth(method="lm", linetype = "dashed", size=0.5, colour="black", alpha=0.3) +
  geom_point(aes(fill=population, shape=type), size=3, alpha=0.8) + 
  labs(x="MDS2", y="Dingo Ancestry") + 
  scale_shape_manual(name="Type", values=ancient_modern_shape) +
  scale_fill_manual(name="Population", values=dingo_ngsd_pops_colours) +
  scale_colour_manual(name="Sample Origin", values=c("New"="red", "Published"="black")) +
  guides(fill="none", shape="none", colour="none") + 
  theme_classic()


world <- ne_coastline(scale="medium")

mds_n_metadata_extra_pop_labels <- left_join(mds_n_metadata_extra_pop, metadata, by=join_by("rn"=="plink_name")) # %>%
  # mutate(sample_id_in_analysis := case_when(
  #   Stagger == "Curracurrang" ~ "Curracurrang.2k",
  #   Stagger == "Nullarbor" ~ "Nullarbor.1k",
  #   .default = sample_id_in_analysis
  # ))
  

new_sample_map <- ggplot() +
  geom_sf(data=world, fill="white") +
  geom_point(data=mds_n_metadata_extra_pop, aes(x=Longitude, y=Latitude, fill=population, shape=type, colour=new, size=new), stroke=1, alpha=0.8) + 
  geom_label_repel(data=mds_n_metadata_extra_pop_labels %>% filter(!is.na(sample_id_in_analysis)), aes(x=Longitude, y=Latitude, label=sample_id_in_analysis),
                   nudge_x=1, nudge_y=-2) + 
  geom_label_repel(data=mds_n_metadata_extra_pop_labels %>% filter(Stagger == "Curracurrang") %>% head(n=1), 
                   aes(x=Longitude, y=Latitude, label="Curracurrang.2k"), nudge_x=10, nudge_y=2) + 
  geom_label_repel(data=mds_n_metadata_extra_pop_labels %>% filter(Stagger == "Nullarbor") %>% head(n=1), 
                   aes(x=Longitude, y=Latitude, label="Nullarbor.1k"), nudge_x=0, nudge_y=-2) + 
  coord_sf(xlim=c(155,110), ylim=c(-45, 0)) +
  scale_shape_manual(name="Type", values=ancient_modern_shape) +
  scale_fill_manual(name="Population", values=dingo_ngsd_pops_colours) +
  scale_colour_manual(name="Sample Origin", values=c("New"="red", "Published"="black")) +
  scale_size_manual(name="Sample Origin", values = c("New"=5, "Published"=4)) + 
  guides(fill="none", shape="none", colour="none", size="none") + 
  #guides(fill=guide_legend(override.aes = list(shape = 21)), colour=guide_legend(override.aes = list(shape = 21))) +
  theme_void() +
  theme(
    legend.position = "bottom"
  )
new_sample_map
ggsave("new_sample_map.png", new_sample_map, dpi=500, height = 10, width=10)
#ggplotly()
final_plot <- new_sample_map | (dingo_ngsd_pca/(admixture_pc1|admixture_pc2)) 
ggsave("figure_dingo_pca.png", final_plot, dpi=500, height=7, width=12)
  
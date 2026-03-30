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
  scales,
  data.table,
  rnaturalearth,
  plotly,
  stringr,
  patchwork,
  ggrain,
  Ternary,
  ggtern
)

setwd("~/Library/CloudStorage/Box-Box/projects/dingo/new_dingo/outgroup_f3")


############ Historical Dog in Australia...




# ################## f4
# data = "merge_masked_all_qpdstat"
# modern_dingoes <- fread("modern_dingoes.list", header=F, col.names = "id")$id
# ancient_dingoes <- fread("ancient_dingoes.list", header=F, col.names = "id")$id
# sahul_dogs <- fread("sahul_dogs.list", header=F, col.names = "id")$id
# ancient_dingoes_f3 <- fread("ancient_dingoes_for_f3.list", header=F, col.names = "id")$id

pca_samples_only_dingo_treetest <- fread("../data/pca_samples_only_dingo_treetest.list", header=F, col.names = c("id"))$id

f4_matrix_treetest <- qpdstat("merge_masked_all_v3_treetest",
                              pop1="CoyoteCalifornia", 
                              pop2="Papua_Indonesia_518bp", 
                              pop3="Dingo_Ancient_Nullarbor",
                              pop4=c(pca_samples_only_dingo_treetest, "HighlandWildDog01", "HighlandWildDog02", "HighlandWildDog03", "NewGuineaSingingDog11", "NewGuineaSingingDog12", "NewGuineaSingingDog15", "NewGuineaSingingDog17"),
                              verbose = TRUE, 
                              sure=TRUE,
                              allsnps=TRUE,
                              auto_only=FALSE,
                              boot=F,
                              adjust_pseudohaploid = TRUE
)
write.csv(f4_matrix_treetest, "f4_matrix_treetest_v3.csv", row.names = F, quote = F)


f4_matrix_treetest_2 <- qpdstat("merge_masked_all_v3_treetest_2",
                                pop1="CoyoteCalifornia", 
                                pop2=pca_samples_only_dingo_treetest, 
                                pop3="HighlandWildDog",
                                pop4="Papua_Indonesia_518bp",
                                verbose = TRUE, 
                                sure=TRUE,
                                allsnps=TRUE,
                                auto_only=FALSE,
                                boot=F,
                                adjust_pseudohaploid = TRUE
)
write.csv(f4_matrix_treetest_2, "f4_matrix_treetest_v3_2.csv", row.names = F, quote = F)

f4_matrix_treetest_3 <- qpdstat("merge_masked_all_v3_treetest_2",
                                pop1="CoyoteCalifornia", 
                                pop2=pca_samples_only_dingo_treetest, 
                                pop3="HighlandWildDog",
                                pop4="NewGuineaSingingDog",
                                verbose = TRUE, 
                                sure=TRUE,
                                allsnps=TRUE,
                                auto_only=FALSE,
                                boot=F,
                                adjust_pseudohaploid = TRUE
)
write.csv(f4_matrix_treetest_3, "f4_matrix_treetest_v3_3.csv", row.names = F, quote = F)

f4_matrix_treetest_4 <- qpdstat("merge_masked_all_v3_treetest_2",
                                pop1="CoyoteCalifornia", 
                                pop2=pca_samples_only_dingo_treetest, 
                                pop3="Papua_Indonesia_518bp",
                                pop4="NewGuineaSingingDog",
                                verbose = TRUE, 
                                sure=TRUE,
                                allsnps=TRUE,
                                auto_only=FALSE,
                                boot=F,
                                adjust_pseudohaploid = TRUE
)
write.csv(f4_matrix_treetest_4, "f4_matrix_treetest_v3_4.csv", row.names = F, quote = F)

#### Checking affinity to Curracurrang
### f4(X, curracurrang; AndeaxFox, Y)
pca_affinity_samples <- fread("../data/pca_samples_dingo_affinity.list", header=F, col.names = c("id"))$id

f4_matrix_curracurrang_affinity <- qpdstat("merge_masked_all_v3", 
                                           pop1="CoyoteCalifornia", 
                                           pop2=pca_affinity_samples, 
                                           pop3="AndeanFox", 
                                           pop4="Dingo_Ancient_Curracurrang",
                                           verbose = TRUE, 
                                           sure=TRUE,
                                           allsnps=TRUE,
                                           auto_only=FALSE,
                                           boot=FALSE,
                                           adjust_pseudohaploid = TRUE,
                                           qpfstats=T)

#### Checking affinity to Nullarbor
f4_matrix_nullarbor_affinity <- qpdstat("merge_masked_all_v3", 
                                        pop1="CoyoteCalifornia", 
                                        pop2=pca_affinity_samples, 
                                        pop3="AndeanFox", 
                                        pop4="Dingo_Ancient_Nullarbor",
                                        verbose = TRUE, 
                                        sure=TRUE,
                                        allsnps=TRUE,
                                        auto_only=FALSE,
                                        boot=FALSE,
                                        adjust_pseudohaploid = T,
                                        qpfstats=T)

pca_affinity_ngsd_samples <- fread("../data/pca_samples_ngsd_affinity.list", header=F, col.names = c("id"))$id
f4_matrix_ngsd_affinity <- qpdstat("merge_masked_all_v3_pap_affinity", 
                                   pop1="CoyoteCalifornia", 
                                   pop2=pca_affinity_ngsd_samples, 
                                   pop3="AndeanFox", 
                                   pop4="Papua_Indonesia_518bp",
                                   verbose = TRUE, 
                                   sure=TRUE,
                                   allsnps=TRUE,
                                   auto_only=FALSE,
                                   boot=FALSE,
                                   adjust_pseudohaploid = TRUE)

write.csv(f4_matrix_curracurrang_affinity, "f4_matrix_curracurrang_affinity_v3.csv", quote=F, row.names = F)
write.csv(f4_matrix_nullarbor_affinity, "f4_matrix_nullarbor_affinity_v3.csv", quote=F, row.names = F)
write.csv(f4_matrix_ngsd_affinity, "f4_matrix_ngsd_affinity_v3.csv", quote=F, row.names = F)

metadata <- fread("../metadata/metadata_for_wgs_samples.tsv") %>%
  janitor::clean_names()

affinity_f4 <- left_join(f4_matrix_curracurrang_affinity %>%
                           dplyr::rename(curracurrang_est = est,
                                         curracurrang_se = se,
                                         curracurrang_z = z,
                                         curracurrang_n=n) %>%
                           dplyr::select(pop2, curracurrang_est, curracurrang_se, curracurrang_z, curracurrang_n),
                         f4_matrix_nullarbor_affinity %>%
                           dplyr::rename(nullarbor_est = est,
                                         nullarbor_se = se,
                                         nullarbor_z = z,
                                         nullarbor_n=n) %>%
                           dplyr::select(pop2, nullarbor_est, nullarbor_se, nullarbor_z, nullarbor_n),
                         by="pop2") %>%
  mutate(affinity = (-nullarbor_est+curracurrang_est)/(curracurrang_est+nullarbor_est))

all_samples <- unique(c(f4_matrix_curracurrang_affinity$pop2, f4_matrix_nullarbor_affinity$pop2, f4_matrix_ngsd_affinity$pop2, f4_matrix_treetest$pop4))


all_sample_coords <- fread("../data/all_sample_coords.tsv")

f4_summary_for_all <- metadata %>%
  filter(sample %in% all_samples,
         !sample %in% c("Foreshore.715bp", "Foreshore.697bp"),
         !str_detect(sample, "NGSD")) %>%
  left_join(., affinity_f4, by=join_by("sample"=="pop2")) %>%
  left_join(., f4_matrix_ngsd_affinity %>% dplyr::rename(ngsd_est = est, ngsd_se=se, ngsd_z=z, ngsd_n=n) %>% dplyr::select(pop2, ngsd_est, ngsd_se, ngsd_z, ngsd_n), by=join_by("sample"=="pop2")) %>%
  left_join(., f4_matrix_treetest %>% filter(pop2 == "Papua_Indonesia_518bp" ) %>% dplyr::rename(tree_est = est, tree_se=se, tree_z=z, tree_n=n) %>% dplyr::select(pop4, tree_est, tree_se, tree_z, tree_n), by=join_by("sample"=="pop4")) %>%
  left_join(., all_sample_coords, by="sample")
write.csv(f4_summary_for_all, "f4_summary_v3.csv", row.names = F, quote = F)

#####################################################################
tree_test_df <- f4_summary_for_all %>%
  filter(!population_label %in% c("", "Captive")) %>%
  mutate(
    population_label2 := case_when(
      population_label == "Ancient NG" ~ "NGSD",
      population_label == "Ancient West" ~ "West", 
      population_label %in% c("Ancient East", "Ancient South") ~ "Ancient Southeast",
      .default = population_label
    ),
    population_label2 = factor(population_label2, levels=c("West", "North", "Central", "Mallee", "Alpine", "East", "K'gari", "HWD", "NGSD", "Ancient Southeast")),
    population_label = factor(population_label, levels=c("Ancient West", "West", "North", "Central", "Mallee", "Alpine", "East", "K'gari", "HWD", "NGSD", "Ancient NG", "Ancient East", "Ancient South"))
    
  ) 

ancient_modern_affinity_plots <- ggplot() +
  geom_violin(data=tree_test_df, aes(x=population_label, y=affinity, fill=population_label), trim = FALSE, scale = "width", colour=NA, alpha=0.5) +
  geom_jitter(data=tree_test_df %>% filter(age=="Modern"), aes(x=population_label, y=affinity, colour=population_label), width=0.05, size=2, alpha=0.7) +
  geom_jitter(data=tree_test_df %>% filter(age=="Ancient"), aes(x=population_label, y=affinity, fill=population_label, label=sample), width=0.05, size=3, shape=24) +
  labs(x="Population", y="Affinity") +
  scale_fill_manual(values=c(pop_colours)) +
  scale_colour_manual(values=c(pop_colours)) +
  guides(fill="none", shape="none", colour="none") +
  coord_cartesian(ylim=c(-0.03,0.03)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle=90)
  )

ancient_modern_affinity_plots


tree_test_plot <- tree_test_df %>%
  filter(!str_detect(sample,"HighlandWildDog"),
         !str_detect(sample,"NewGuineaSingingDog")) %>%
  arrange(tree_est) %>%
  mutate(sample = factor(sample, levels=sample)) %>%
  filter(tree_est != "") %>%
  ggplot(aes(x=sample, y=tree_est)) +
  geom_hline(yintercept = 0, alpha=0.1) +
  geom_errorbar(aes(ymin=tree_est-tree_se, ymax=tree_est+tree_se, colour=abs(tree_z)), width=0.5) +
  geom_point(aes(x=sample, y=tree_est, fill=population_label, shape=age), size=4) +
  labs(x="X", y=expression(italic(f[4](CoyoteCalifornia~","~Ancient~New~Guinean~Dingo~";"~Ancient~Nullarbor~Dingo~","~X)))) +
  #scale_colour_manual(values=c("red", "green")) +
  scale_colour_gradientn(name= "|Z| score",colours=c("green","yellow"), values=rescale(x=c(0, 3), from=c(0,max(tree_test_df$tree_z, na.rm=T))), breaks=c(0.1, 3, 13),labels = c("0","< 3", "14"), na.value = "red") +
  scale_fill_manual(values=pop_colours) +
  scale_shape_manual(values=c("Modern"=21, "Ancient"=24)) + 
  guides(fill="none", shape="none") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

ggplotly(tree_test_df %>%
           filter(!str_detect(sample,"HighlandWildDog"),
                  !str_detect(sample,"NewGuineaSingingDog")) %>%
           arrange(tree_est) %>%
           mutate(sample = factor(sample, levels=sample)) %>%
           filter(tree_est != "") %>%
           ggplot(aes(x=sample, y=tree_est)) +
           geom_point(aes(x=sample, y=tree_est, fill=population_label, shape=age, label=sample), size=4))

tree_test_plot_2 <- tree_test_df %>%
  filter(!str_detect(sample,"HighlandWildDog"),
         !str_detect(sample,"NewGuineaSingingDog")) %>%
  arrange(tree_est) %>%
  mutate(sample = factor(sample, levels=sample)) %>%
  filter(tree_est != "") %>%
  ggplot(aes(y=sample, x=tree_est)) +
  geom_vline(xintercept = 0, alpha=0.1) +
  geom_errorbarh(aes(xmin=tree_est-tree_se, xmax=tree_est+tree_se, colour=abs(tree_z)), height=0.5) +
  geom_point(aes(y=sample, x=tree_est, fill=population_label, shape=age), size=4) +
  labs(y="X", x=expression(italic(f[4](CoyoteCalifornia~","~Ancient~New~Guinean~Dingo~";"~Ancient~Nullarbor~Dingo~","~X)))) +
  #scale_colour_manual(values=c("red", "green")) +
  scale_colour_gradientn(name= "|Z| score",colours=c("green","yellow"), limit=c(0,3), na.value = "red") +
  scale_fill_manual(values=pop_colours) +
  scale_shape_manual(values=c("Modern"=21, "Ancient"=24)) + 
  guides(fill="none", shape="none") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

tree_test_plot_2

ggsave("f4_tree_test.png", tree_test_plot_2, dpi=300, width=8, height=9)

ternary_plot <- ggplot(f4_summary_for_all,
       aes(
         x = nullarbor_est,
         y = curracurrang_est,
         z = ngsd_est,
         shape = age, 
         fill = population_label
       )) +
  #stat_density_tern(aes(colour = population_label ), alpha = .2) +
  geom_point(data=f4_summary_for_all %>% filter(age=="Modern"), aes(group = population_label), size = 3, alpha=0.5) +
  geom_point(data=f4_summary_for_all %>% filter(age=="Ancient"), aes(group = population_label), size = 4, alpha=0.7) +
  scale_shape_manual(values = c(24, 21)) +
  coord_tern() +
  limit_tern(T=.37,L=.37,R=.34) +
  Tlab(label = "", labelarrow ="Ancient Curracurrang\nDingo\n") +
  Llab(label = "", labelarrow = "Ancient Nullarbor\nDingo\n") +
  Rlab(label = "", labelarrow = "\nAncient New Guinean\nDingo") +
  labs(
    shape = "Populations",
    fill = "Longitude",
    colour = "Longitude",
    alpha = NULL
  ) +
  guides(
    alpha = "none",
    color = "none",
    fill="none",
    shape="none") +
  scale_fill_manual(values=pop_colours) + 
  scale_colour_manual(values=pop_colours) + 
  theme_custom(
    base_size = 12,
    tern.plot.background = NULL,
    tern.panel.background = "white",
    col.T = "#FFAE43",
    col.L = "#0D4162",
    col.R = "black",
    col.grid.minor = "white"
  )
ggsave("ternary_plot.png", ternary_plot, dpi=300, width=3, height=3)


f4_affinity_plots <- (ternary_plot | tree_test_plot) + plot_layout(widths = c(1,1))


ggsave("f4_plots.png", f4_affinity_plots, dpi=300, height=6.5, width=12)


tree_test_plot <- ggplot() +
  geom_hline(yintercept = 0, alpha=0.1) +
  geom_violin(data=tree_test_df %>% filter(!population_label %in% c("HWD", "NGSD", "Ancient NG")), aes(x=population_label, y=tree_est), trim = FALSE, scale = "width", alpha=0.5) +
  #geom_errorbar(data=tree_test_df %>% filter(!population_label %in% c("HWD", "NGSD", "Ancient NG")), aes(x=population_label, y=tree_est, ymin=tree_est-tree_se, ymax=tree_est+tree_se), width=0.05, size=2, alpha=0.2, position=position_dodge(width=0.05)) +
  geom_pointrange(data=tree_test_df %>% filter(!population_label %in% c("HWD", "NGSD", "Ancient NG"), age=="Modern"), aes(x=population_label, y=tree_est, colour=as.factor(tree_z>3)), width=0.05, size=2, alpha=0.7, position=position_dodge(width=0.05)) +
  #geom_point(data=tree_test_df %>% filter(!population_label %in% c("HWD", "NGSD", "Ancient NG"), age=="Ancient"), aes(x=population_label, y=tree_est, fill=as.factor(tree_z>3), label=sample), width=0.05, size=3, shape=24,position=position_dodge(width=0.05)) +
  labs(x="Population", y=expression(italic(f[4](CoyoteCalifornia~","~Ancient~New~Guinean~Dingo~";"~AndeanFox~","~X)))) +
  scale_fill_manual(values=c("red", "green")) +
  scale_colour_manual(values=c("red", "green")) +
  guides(fill="none", shape="none", colour="none") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle=90)
  )
tree_test_plot

##### Maybe good to try with longitude...
ggplotly(f4_summary_for_all %>%
  arrange(tree_est) %>%
  mutate(sample=factor(sample, levels=sample)) %>%
  ggplot(aes(x=sample, y=tree_est, shape=as.factor(abs(tree_z)>3), colour=population_label)) +
  geom_errorbar(aes(ymin=tree_est-tree_se, ymax=tree_est+tree_se)) +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)))


f4_summary_for_all %>%
  filter(!is.na(tree_est))


f4_matrix_ancient_dog <- fread("f4_matrix_ancient_dog.csv") %>%
  group_by(pop4) %>%
  arrange(desc(est)) %>%
  dplyr::slice_head(n=20)
  
fread("f4_matrix_ancient_dog.csv") %>%
  filter(pop4 == "Meredith.99bp",
         pop2 != "AfricanHuntingDogKenya") %>%
  arrange(est) %>%
  mutate(pop2 = factor(pop2, levels=pop2)) %>%
  ggplot(aes(x=est, y=pop2)) +
  geom_errorbar(aes(xmin=est-se, xmax=est+se), width=0.2) + 
  geom_point() 


dog1 <- f4_matrix_ancient_dog %>%
  filter(pop4 == "Castlemaine.99bp") %>%
  arrange(est) %>%
  mutate(pop2 = factor(pop2, levels=pop2)) %>%
  ggplot(aes(x=est, y=pop2)) +
  geom_errorbar(aes(xmin=est-se, xmax=est+se), width=0.2) + 
  geom_point() +
  labs(x=expression(italic(f[4]-estimate)), y="Dog Breed") + 
  facet_grid(~pop4) +
  theme_bw()+
  theme(panel.grid = element_blank())

dog2 <- f4_matrix_ancient_dog %>%
  filter(pop4 == "FernCave.100bp") %>%
  arrange(est) %>%
  mutate(pop2 = factor(pop2, levels=pop2)) %>%
  ggplot(aes(x=est, y=pop2)) +
  geom_errorbar(aes(xmin=est-se, xmax=est+se), width=0.2) + 
  geom_point() +
  labs(x=expression(italic(f[4]-estimate)), y="Dog Breed") + 
  facet_grid(~pop4) +
  theme_bw()+
  theme(panel.grid = element_blank())

dog3 <- f4_matrix_ancient_dog %>%
  filter(pop4 == "Meredith.99bp") %>%
  arrange(est) %>%
  mutate(pop2 = factor(pop2, levels=pop2)) %>%
  ggplot(aes(x=est, y=pop2)) +
  geom_errorbar(aes(xmin=est-se, xmax=est+se), width=0.2) + 
  geom_point() +
  labs(x=expression(italic(f[4]-estimate)), y="Dog Breed") + 
  facet_grid(~pop4) +
  theme_bw()+
  theme(panel.grid = element_blank())
dog1 | dog2 | dog3

########## f4_matrix_modern_dog_bits.csv
f4_matrix_modern_dog <- fread("../data/f4_matrix_modern_dog_bits.csv") 
f4_matrix_ancient_dog <- fread("../data/f4_matrix_ancient_dog_bits.csv") %>%
  filter(pop2 != "FernCave.100bp", n > 10000)
dog_cata <- fread("../data/dogs.list")


f4_matrix_ancient_dog %>%
  filter(pop2 != "FernCave.100bp") %>%
  ggplot(aes(x=n, y=after_stat(count), fill=pop2))  +
  geom_histogram()

f4_matrix_modern_dog %>%
  filter(n > 50000) %>%
  rbind(., f4_matrix_ancient_dog) %>%
  left_join(., metadata, by=join_by("pop4"=="sample")) %>%
  left_join(., dog_cata, by=join_by("pop2"=="Breed")) %>%
  filter(population_label == "Alpine") %>%
  group_by(pop2) %>%
  arrange(desc(est)) %>%
  dplyr::slice_head(n=20) %>%
  mutate(pop2 = factor(pop2, levels=unique(pop2))) %>%
  ggplot(aes(x=pop2, y=est)) +
  geom_point() +
  geom_errorbar(aes(xmin=est-se, xmax=est+se), width=0.2) +
  facet_grid(rows = vars(pop4)) +
  theme(
    axis.text.x = element_text(angle=90)
  )
  

f4_matrix_modern_dog_filt <- f4_matrix_modern_dog %>%
  filter(n > 50000) %>%
  group_by(pop4) %>%
  arrange(desc(est)) %>%
  dplyr::slice_head(n=1) %>%
  left_join(., metadata, by=join_by("pop4"=="sample")) %>%
  left_join(., dog_cata, by=join_by("pop2"=="Breed")) %>%
  mutate(Origin := case_when(
    is.na(Origin) ~ "North America",
    .default = Origin
  )) %>%
  ungroup()

class_table <- table(f4_matrix_modern_dog_filt %>%
  dplyr::select(population_label, Class) %>%
  group_by(population_label))

class_table <- f4_matrix_modern_dog_filt %>%
  filter(population_label != "") %>%
  count(population_label, Class)

origin_table <- f4_matrix_modern_dog_filt %>%
  filter(population_label != "") %>%
  count(population_label, Origin)

breed_plot <- ggplot(class_table, aes(x = population_label, y = n, fill = Class)) +
  geom_col() +
  scale_fill_manual(values=c("#264653","#287271","#2a9d8f","#8ab17d","#e9c46a","#efb366","#f4a261","#ee8959","#e76f51")) +
  theme_bw() +
  labs(x = "Population",
       y = "Count",
       fill = "Class") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        legend.position = "bottom")

origin_plot <- ggplot(origin_table, aes(x = population_label, y = n, fill = Origin)) +
  geom_col() +
  theme_bw() +
  scale_fill_manual(values=c("#006954","#0e5e7a","#f7915e","#ffd374","#063447","#963831","#4a3941","#b9d29d")) + 
  labs(x = "Population",
       y = "Count",
       fill = "Origin") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        legend.position = "bottom")

ancient_dogs <- ggarrange(dog1, dog2 , dog3, ncol=3, labels = c("a", "b", "c"))

modern_dogs <- ggarrange(breed_plot, origin_plot, ncol=2, legend = "bottom", labels=c("d", "e"))

dog_summary <- ggarrange(ancient_dogs, modern_dogs, ncol=1)
ggsave("f4_dog_summary.png", dog_summary, dpi=300, height=10, width=15)


ggplot() +
  geom_sf(data=australia)



##################
# 
# final_sample_list <- fread("~/Library/CloudStorage/Box-Box/projects/dingo/new_dingo/samples_to_keep_pca.list",
#                            header = F, col.names = c("id"))$id
# 
# 
# 
# f4_matrix_ngsd_affinity %>%
#   filter(pop2 %in% c(final_sample_list, "HighlandWildDog02", "HighlandWildDog03")) %>%
#   arrange(est) %>%
#   mutate(pop2=factor(pop2, levels=pop2)) %>%
#   ggplot(aes(x=est, y=pop2)) +
#   geom_errorbarh(aes(xmin=est - se, xmax=est+se)) +
#   geom_point() +
#   theme_bw()
# 
# affinity_f4 <- left_join(f4_matrix_curracurrang_affinity %>%
#                            dplyr::rename(curracurrang_est = est,
#                                          curracurrang_se = se,
#                                          curracurrang_z = z) %>%
#                            dplyr::select(pop2, curracurrang_est, curracurrang_se, curracurrang_z),
#                          f4_matrix_nullarbor_affinity %>%
#                            dplyr::rename(nullarbor_est = est,
#                                          nullarbor_se = se,
#                                          nullarbor_z = z) %>%
#                            dplyr::select(pop2, nullarbor_est, nullarbor_se, nullarbor_z),
#                          by="pop2") %>%
#   mutate(affinity = (-nullarbor_est+curracurrang_est)/(curracurrang_est+nullarbor_est))
# 
# write.csv(affinity_f4, "nulla_curra_affinity_f4_new.csv", quote = F, row.names = F)
# 
# ggplotly(affinity_f4 %>%
#   filter(pop2 %in% c(final_sample_list, "HighlandWildDog02", "HighlandWildDog03")) %>%
#   arrange(affinity) %>%
#   mutate(pop2=factor(pop2, levels=pop2)) %>%
#   ggplot(aes(x=affinity, y=pop2)) +
#   geom_point() +
#   theme_bw())
# 
# ### Remove related samples 
# related_samples <- fread("../king/wgs_samples_to_remove.list", header = F, col.names = c("id"))$id
# 
# 
# ggplotly(affinity_f4 %>%
#   filter(!pop2 %in% c(related_samples, "Foreshore.697bp", "Foreshore.715bp")) %>%
#   ggplot(aes(x=curracurrang_est, y=nullarbor_est, label=pop2)) +
#   geom_abline(intercept = 0, alpha=0.7, linetype="dashed", size=0.5) + 
#   geom_errorbarh(aes(xmin=curracurrang_est - curracurrang_se, xmax=curracurrang_est+curracurrang_se)) + 
#   geom_errorbar(aes(ymin=nullarbor_est - nullarbor_se, ymax=nullarbor_est+nullarbor_se)) + 
#   geom_point(aes(label=pop2)))
# 
# 
# # affinity_plot <- affinity_f4 %>%
# #   filter(!pop2 %in% related_samples) %>%
# #   ggplot(aes(x=curracurrang_est, y=nullarbor_est, label=pop2)) +
# #   geom_abline(intercept = 0, alpha=0.7, linetype="dashed", size=0.5) + 
# #   geom_errorbarh(aes(xmin=curracurrang_est - curracurrang_se, xmax=curracurrang_est+curracurrang_se)) + 
# #   geom_errorbar(aes(ymin=nullarbor_est - nullarbor_se, ymax=nullarbor_est+nullarbor_se)) + 
# #   geom_point() +
# #   scale_fill_manual(name="Population", values=dingo_ngsd_pops_colours) +
# #   scale_colour_manual(name="Sample Origin", values=c("New"="red", "Published"="black")) +
# #   scale_size_manual(name="Sample Origin", values = c("New"=4, "Published"=3)) + 
# #   scale_shape_manual(name="Type", values=ancient_modern_shape) +
# #   scale_alpha_manual(name="Type", values=c("Ancient"=1, "Modern"=0.8)) + 
# #   labs(x="f4(CoyoteCalifornia, Curracurrang.2k; AndeanFox, X)", y="f4(CoyoteCalifornia, Nullarbor.1k; AndeanFox, X)") + 
# #   guides(fill="none", colour="none",
# #          size="none", alpha="none", shape="none") +
# #   theme_classic()
# 
# 
# #### Affinity of modern dingo to ancients...
# modern_dingoes <- pops_f4_dingoes %>%
#   filter(!str_detect(pop, "Dingo_Ancient"))
# 
# ancient_dingoes <- pops_f4_dingoes %>%
#   filter(str_detect(pop, "Dingo_Ancient"),
#          !pop %in% c("Dingo_Ancient_Gatton_LP_14", "Dingo_Ancient_VIC_04", "Dingo_Ancient_VIC_05", "Dingo_Ancient_VIC_10", "Dingo_Ancient_VIC_LP_01_02.2807bp"))
# 
# f4_matrix_modern_ancient_affinity <- qpdstat(data, 
#                                              pop1="CoyoteCalifornia", 
#                                              pop2=c(ancient_dingoes$pop, 
#                                                     "Papua_Indonesia_518bp", 
#                                                     "Dingo_Ancient_Nullarbor", 
#                                                     "Dingo_Ancient_Curracurrang"), 
#                                              pop3="AndeanFox", 
#                                              pop4=modern_dingoes$pop,
#                                              verbose = TRUE, 
#                                              sure=TRUE,
#                                              allsnps=TRUE,
#                                              auto_only=FALSE,
#                                              boot=FALSE, f4mode=TRUE)
# 
# write.csv(f4_matrix_modern_ancient_affinity, "modern_to_ancient_affinity_f4.csv", quote = F, row.names = F)
# 
# 
# modern_dingo_affinity <- f4_matrix_modern_ancient_affinity %>%
#   filter(!pop4 %in% c("NewGuineaSingingDog_BNGS2", "NewGuineaSingingDog_BNGS1", "Dingo_Alpine_722g")) %>%
#   left_join(., mds_n_metadata_extra_pop, by=join_by("pop4"=="rn")) %>%
#   left_join(., metadata, by=join_by("pop2"=="plink_name")) %>%
#   left_join(., metadata_old_samples, by=join_by("pop2"=="plink_name")) %>%
#   mutate(sample_id_in_analysis := case_when(
#     str_detect(pop2, "Curracurrang") ~ "Curracurrang.2k",
#     str_detect(pop2, "Nullarbor") ~ "Nullarbor.1k",
#     .default = sample_id_in_analysis
#   )) %>%
#   group_by(sample_id_in_analysis) %>%
#   arrange(sample_id_in_analysis, est) %>%
#   ggplot(aes(x=est, y=pop4)) +
#   geom_errorbarh(aes(xmin=est-se, xmax=est+se), height = .1) + 
#   geom_point(aes(shape=type, fill=population), size=3) +
#   facet_grid(~sample_id_in_analysis) +
#   scale_fill_manual(name="Population", values=dingo_ngsd_pops_colours) +
#   scale_shape_manual(name="Type", values=ancient_modern_shape) +
#   labs(x="f4(CoyoteCalifornia, Y; AndeanFox, X)") + 
#   guides(fill = guide_legend(override.aes = list(shape = 21, size=3)), shape="none") + 
#   theme_bw() +
#   theme(
#     axis.text.y = element_blank(),
#     axis.title.y = element_blank(),
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#     legend.position = "bottom"
#   )
# 
# ggsave("modern_dingo_affinity.png", modern_dingo_affinity, dpi=500, height = 5, width=16)
# 
# 
# ### Checking NGSD affinity
# f4_matrix_ngsd_affinity <- qpdstat(data, pop1="CoyoteCalifornia", pop2="Papua_Indonesia_518bp", pop3="AndeanFox", pop4=pops_f4,
#                           verbose = TRUE, 
#                           sure=TRUE,
#                           allsnps=TRUE,
#                           auto_only=FALSE,
#                           boot=FALSE,
#                           f4mode=TRUE)
# 
# 
# f4_tests <- rbind(f4_matrix_curracurrang_affinity, )
# 
# metadata_old_samples <- fread("../metadata/metadata_old_samples.tsv") %>%
#   janitor::clean_names()
# 
# ngsd_affinity_plot <- f4_matrix_ngsd_affinity %>%
#   filter(!pop4 %in% c("NewGuineaSingingDog_BNGS2", "NewGuineaSingingDog_BNGS1", "Dingo_Ancient_Gatton_LP_14",
#                       "Dingo_Ancient_VIC_04", "Dingo_Ancient_VIC_05", "Dingo_Ancient_VIC_10",
#                       "Dingo_Ancient_Foreshore_697bp", "Dingo_Ancient_Foreshore_715bp", "Dingo_Alpine_722g",
#                       "Dingo_Ancient_VIC_LP_01_02.2807bp")) %>%
#   left_join(., mds_n_metadata_extra_pop, by=join_by("pop4"=="rn")) %>%
#   left_join(., metadata, by=join_by("pop4"=="plink_name")) %>%
#   left_join(., metadata_old_samples, by=join_by("pop4"=="plink_name")) %>%
#   mutate(sample_id_in_analysis := case_when(
#     str_detect(pop4, "Curracurrang") ~ paste0("Curracurrang.", date.y, "bp"),
#     str_detect(pop4, "Nullarbor") ~ paste0("Nullarbor.", date.y, "bp"),
#     str_detect(pop4, "FraserIsland") ~ str_replace(pop4, "Dingo_FraserIsland", "Kgari"),
#     is.na(sample_id_in_analysis) & str_detect(pop4, "Dingo_") ~ str_replace(pop4, "Dingo_", ""),
#     str_detect(pop4, "NewGuinea") ~ pop4,
#     .default = sample_id_in_analysis
#   )) %>%
#   arrange(desc(est)) %>% 
#   mutate(sample_id_in_analysis = factor(sample_id_in_analysis, levels=sample_id_in_analysis)) %>%
#   ggplot(aes(y=sample_id_in_analysis, x=est, label=pop4, fill=population, shape=type, alpha=type)) +
#   geom_errorbarh(aes(xmin=est - se, xmax=est+se)) + 
#   geom_point(aes(size=new), size=4) +
#   scale_fill_manual(name="Population", values=dingo_ngsd_pops_colours) +
#   scale_colour_manual(name="Sample Origin", values=c("New"="red", "Published"="black")) +
#   scale_size_manual(name="Sample Origin", values = c("New"=4, "Published"=3)) + 
#   scale_shape_manual(name="Type", values=ancient_modern_shape) +
#   scale_alpha_manual(name="Type", values=c("Ancient"=1, "Modern"=0.8)) + 
#   labs(y="Sample Name", x="f4(CoyoteCalifornia, Papua.518bp; AndeanFox, X)") + 
#   guides(fill=guide_legend(override.aes = list(shape = 21, size=5)), colour="none",
#          size="none", alpha="none", shape="none") +
#   theme_classic()
# 
# ngsd_affinity_plot_2 <- f4_matrix_ngsd_affinity %>%
#   filter(!pop4 %in% c("NewGuineaSingingDog_BNGS2", "NewGuineaSingingDog_BNGS1", "Dingo_Ancient_Gatton_LP_14",
#                       "Dingo_Ancient_VIC_04", "Dingo_Ancient_VIC_05", "Dingo_Ancient_VIC_10",
#                       "Dingo_Ancient_Foreshore_697bp", "Dingo_Ancient_Foreshore_715bp", "Dingo_Alpine_722g",
#                       "Dingo_Ancient_VIC_LP_01_02.2807bp")) %>%
#   left_join(., mds_n_metadata_extra_pop, by=join_by("pop4"=="rn")) %>%
#   left_join(., metadata, by=join_by("pop4"=="plink_name")) %>%
#   left_join(., metadata_old_samples, by=join_by("pop4"=="plink_name")) %>%
#   mutate(sample_id_in_analysis := case_when(
#     str_detect(pop4, "Curracurrang") ~ paste0("Curracurrang.", date.y, "bp"),
#     str_detect(pop4, "Nullarbor") ~ paste0("Nullarbor.", date.y, "bp"),
#     str_detect(pop4, "FraserIsland") ~ str_replace(pop4, "Dingo_FraserIsland", "Kgari"),
#     is.na(sample_id_in_analysis) & str_detect(pop4, "Dingo_") ~ str_replace(pop4, "Dingo_", ""),
#     str_detect(pop4, "NewGuinea") ~ pop4,
#     .default = sample_id_in_analysis
#   ),
#   date := case_when(
#     !is.na(date.x) ~ date.x,
#     !is.na(date.y) ~ date.y,
#     sample_id_in_analysis == "Gatton.nd" ~ 2500,
#     .default = 0
#   )) %>%
#   rowwise() %>%
#   mutate(
#     date := case_when(
#       date == 0 ~ sample(c(0,-20,-40,-60,-80,-100,-120,-140,-160,-180,-200),1),
#       .default = date
#     )
#   ) %>%
#   arrange(desc(est)) %>% 
#   mutate(sample_id_in_analysis = factor(sample_id_in_analysis, levels=sample_id_in_analysis)) %>%
#   ggplot(aes(x=date, y=est, label=pop4, fill=population, shape=type)) +
#   geom_rect(aes(xmin=-250,xmax=50, ymin=-Inf, ymax=Inf), fill="grey90", alpha=0.1) + 
#   geom_rect(aes(xmin=2450,xmax=2550, ymin=-Inf, ymax=Inf), fill="grey90", alpha=0.1) + 
#   geom_errorbar(aes(ymin=est - se, ymax=est+se)) + 
#   geom_point(aes(size=new), size=4) +
#   scale_fill_manual(name="Population", values=dingo_ngsd_pops_colours) +
#   scale_colour_manual(name="Sample Origin", values=c("New"="red", "Published"="black")) +
#   scale_size_manual(name="Sample Origin", values = c("New"=4, "Published"=3)) + 
#   scale_shape_manual(name="Type", values=ancient_modern_shape) +
#   labs(y="f4(CoyoteCalifornia, Papua.518bp; AndeanFox, X)", x="Date (BP)") + 
#   guides(fill=guide_legend(override.aes = list(shape = 21, size=5)), colour="none",
#          size="none", alpha="none", shape="none") +
#   scale_x_continuous(breaks=c(-100, 500, 1000, 1500, 2000, 2500), labels = c("0", "500", "1000", "1500", "2000", "ND")) +
#   theme_classic()
# ngsd_affinity_plot_2
# combined_plot <- (affinity_plot | ngsd_affinity_plot) / ngsd_affinity_plot_2
# ggsave("combined_affinities.png", combined_plot, dpi=500, height=5, width=13)
# 
# 
# write.csv(f4_matrix_CPAX, "f4_pops_Coyote_Papua_Andean_x.csv", quote = F, row.names = F)
# 
# all_f4_tests <- rbind(f4_matrix_curracurrang_affinity, f4_matrix_nullarbor_affinity) %>%
#   rbind(., f4_matrix_CPAX)
# 
# write.csv(all_f4_tests, "all_f4_tests.csv", quote = F, row.names = F)
# 
# 
# 
# ##### Admixture F4 
# f4_matrix_CPAncientX <- qpdstat(data, pop1="CoyoteCalifornia", pop2="Papua_Indonesia_518bp", pop3=c("Zhokhov9500BP.CGG6", "Germany_7k", "Ireland_Neolithic.Newgrange", "Russia_Baikal_7k"), pop4=pops_f4,
#                                 verbose = TRUE, 
#                                 sure=TRUE,
#                                 allsnps=TRUE,
#                                 auto_only=FALSE,
#                                 boot=FALSE,
#                                 f4mode=TRUE)
# 
# write.csv(f4_matrix_CPAncientX, "f4_pops_Coyote_Papua_Ancient_Dogs_x.csv", quote = F, row.names = F)
# f4_matrix_CPAncientX <- fread("f4_pops_Coyote_Papua_Ancient_Dogs_x.csv")
# 
# gene_flow_ngsd <- f4_matrix_CPAncientX %>%
#   filter(!pop4 %in% c("NewGuineaSingingDog_BNGS2", "NewGuineaSingingDog_BNGS1", "Dingo_Ancient_Gatton_LP_14",
#                       "Dingo_Ancient_VIC_04", "Dingo_Ancient_VIC_05", "Dingo_Ancient_VIC_10",
#                       "Dingo_Ancient_Foreshore_697bp", "Dingo_Ancient_Foreshore_715bp", "Dingo_Alpine_722g")) %>%
#   left_join(., mds_n_metadata_extra_pop, by=join_by("pop4"=="rn")) %>%
#   arrange(desc(est)) %>% 
#   mutate(pop4 = factor(pop4, levels=unique(pop4))) %>%
#   ggplot(aes(y=pop4, x=est, fill=population, shape=type)) +
#   geom_errorbarh(aes(xmin=est - se, xmax=est+se), height=.1) + 
#   geom_point(size=3) +
#   facet_grid(~pop3) + 
#   scale_fill_manual(name="Population", values=dingo_ngsd_pops_colours) +
#   scale_colour_manual(name="Sample Origin", values=c("New"="red", "Published"="black")) +
#   scale_shape_manual(name="Type", values=ancient_modern_shape) +
#   labs(y="X", x="f4(CoyoteCalifornia, NGSD_518; Y, X)") + 
#   guides(fill="none", colour="none",
#          size="none", alpha="none", shape="none") +
#   theme_classic()
# 
# ggsave("gene_flow_ngsd.png", gene_flow_ngsd, dpi=500, height=5, width=12)
# 
# 
# 
# f4_matrix_CBGX <- qpdstat(data, pop1="CoyoteCalifornia", pop2="Russia_Baikal_7k", pop3="GermanShepherdDog", pop4=pops_f4,
#                           verbose = TRUE, 
#                           sure=TRUE,
#                           allsnps=TRUE,
#                           auto_only=FALSE,
#                           boot=FALSE,
#                           f4mode=TRUE)
# 
# write.csv(f4_matrix_CBGX, "f4_pops_Coyote_Baikal_German_x.csv", quote = F, row.names = F)
# 
# f4_matrix_CBWX <- qpdstat(data, pop1="CoyoteCalifornia", pop2="WolfSyria", pop3="GermanShepherdDog", pop4=pops_f4,
#                           verbose = TRUE, 
#                           sure=TRUE,
#                           allsnps=TRUE,
#                           auto_only=FALSE,
#                           boot=FALSE,
#                           f4mode=TRUE)
# write.csv(f4_matrix_CBWX, "f4_pops_Coyote_Wolf_German_x.csv", quote = F, row.names = F)
# 
# 
# ###### NGSD plot
# bad_samples_2 <- c("NewGuineaSingingDog_BNGS1", "NewGuineaSingingDog_BNGS2")
# f4_plot_CPGX <- f4_matrix_CPGX %>% 
#   filter(!pop4 %in% bad_samples_2) %>%
#   arrange(est) %>%
#   mutate(pop4 = factor(pop4, levels=pop4)) %>%
#   ggplot(aes(x=pop4, y=est)) +
#   geom_point(shape=21, size=2) +
#   geom_errorbar(aes(ymin=est-se, ymax=est+se)) +
#   labs(x="f4(CoyoteCalifornia,NGSD_518bp;GermanShepherdDog,X)")+
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#     panel.grid.major.x  = element_blank(),
#   )
# 
# ggsave("f4_plot_CPGX.png", f4_plot_CPGX, dpi=500)
# 
# f4_plot_CPAncientX <- f4_matrix_CPAncientX %>% 
#   filter(!pop4 %in% bad_samples_2) %>%
#   arrange(est) %>%
#   mutate(pop4 = factor(pop4, levels=unique(pop4))) %>%
#   ggplot(aes(x=pop4, y=est)) +
#   geom_point(shape=21, size=2) +
#   facet_grid(~pop3) + 
#   geom_errorbar(aes(ymin=est-se, ymax=est+se)) +
#   labs(x="f4(CoyoteCalifornia,NGSD_518bp;Y,X)")+
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#     panel.grid.major.x  = element_blank(),
#   )
# f4_plot_CPAncientX
# 
# ggsave("f4_plot_CPAncientX.png", f4_plot_CPAncientX, dpi=500, width=15, height = 5)
# 
# f4_matrix_CBWX %>% 
#   filter(!pop4 %in% bad_samples_2) %>%
#   arrange(est) %>%
#   mutate(pop4 = factor(pop4, levels=pop4)) %>%
#   ggplot(aes(x=pop4, y=est)) +
#   geom_point(shape=21, size=2) +
#   geom_errorbar(aes(ymin=est-se, ymax=est+se)) +
#   labs(x="f4(CoyoteCalifornia,NGSD_518bp;SyrianWolf,X)")+
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#     panel.grid.major.x  = element_blank(),
#   )
# 
# 
# f4_matrix_CGattonGX <- fread("f4_matrix_CGattonGX.csv")
# 
# # Armenia_Hovk1.HOV4 is a dhole
# # Quebec is very close to our sample for some reason. Removing it 
# 
# remove_outliers <- c("EthiopianWolf", "QuebecWolf")
# 
# 
# f4_plot_CGattonGX <- f4_matrix_CGattonGX %>%
#   filter(!pop4 %in% remove_outliers) %>%
#   arrange(est) %>%
#   mutate(pop4 = factor(pop4, levels=pop4)) %>%
#   ggplot(aes(x=pop4, y=est)) +
#   geom_point(shape=21, size=2) +
#   geom_errorbar(aes(ymin=est-se, ymax=est+se)) +
#   labs(x="f4(CoyoteCalifornia,Gatton;GermanShepherdDog,X)")+
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4),
#     panel.grid.major.x  = element_blank(),
#   )
# 
# ggsave("f4_plot_CGattonGX.png", f4_plot_CGattonGX, dpi=500, height=10, width=15)  
# 
# 
# f4_plot_CXAncient_DholeGatton <- fread("f4_matrix_CXAncient_DholeGatton.csv") %>%
#   arrange(est) %>%
#   mutate(pop2 = factor(pop2, levels=pop2)) %>%
#   ggplot(aes(x=est,y=pop2)) +
#   geom_point(shape=21) +
#   geom_errorbar(aes(colour=as.factor(abs(z) > 3), xmin=est-se, xmax=est+se)) +
#   theme(
#     axis.text.y = element_text(size=4),
#     panel.grid.major.x  = element_blank(),
#   )
# 
# ggsave("f4_plot_CXAncient_DholeGatton.png", f4_plot_CXAncient_DholeGatton, dpi=500, height=20, width=7)
# 
# fread("ff4_matrix_CAncient_DholeGattonX.csv") %>% 
#   arrange(est)  %>%
#   mutate(pop4 = factor(pop4, levels=pop4)) %>%
#   ggplot(aes(x=est,y=pop4)) +
#   geom_point(shape=21) +
#   geom_errorbar(aes(xmin=est-se, xmax=est+se))
# 
# 
# ####### Affinity f4
# # f4(X, CoyoteCalifornia; Y, AndeanFox)
# dingo_pops <- c("")
# ancient_ancestry <- c("")
# ancient_dingo_affinity <- qpdstat(data, pop1=dingo_pops, pop2="CoyoteCalifornia", pop3=ancient_ancestry, pop4="AndeanFox",
#                                   verbose = TRUE, 
#                                   sure=TRUE,
#                                   allsnps=TRUE,
#                                   auto_only=FALSE,
#                                   boot=FALSE,
#                                   f4mode=TRUE)
# 
# 
# 
# ###### What type of dog are the ancient dogs
# f4_matrix_ancient_dog <- fread("f4_matrix_ancient_dog.csv")
# 
# # f4(coyote,X;dhole,gatton)
# # f4(coyote,dhole;Gatton,X)
# 
# # f4(AdeanFox,X;dhole,gatton)
# # f4(AdeanFox,dhole;Gatton,X)
# 
# 
# 
# 
# left_right <- c("Dingo_Ancient_Nullarbor", "GermanShepherdDog", "HighlandWildDog", "Papua_Indonesia_518bp")
# right_fix <- c("CoyoteCalifornia", "Russia_Baikal_7k", "Germany_7k", "Ireland_Neolithic.Newgrange") # , 
# 
# all_combinations_all <- unlist(
#   lapply(seq_len(3), function(x) combn(left_right, x, simplify = FALSE)), 
#   recursive = FALSE
# )
#

#
# all_combinations <- 
#   all_combinations_all[
#     sapply(all_combinations_all, function(x) 
#       "Dingo_Ancient_Nullarbor" %in% x
#     )
#   ][1:6]
# 
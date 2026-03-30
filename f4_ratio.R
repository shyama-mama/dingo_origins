pacman::p_load(
  admixtools,
  dplyr,
  data.table
)


setwd("~/Library/CloudStorage/Box-Box/projects/dingo/new_dingo/f4_ratio")

# f4 ratio
# f4(CoyoteCalifornia,Nullarbor_1k;X,Papua)           - f4(1,2;3,4)
# ---------------------------------
# f4(CoyoteCalifornia,Nullarbor_1k;Curracurrang_2k,Papua)       - f4(1,2;5,4)
# a is between 3 and 5 where 3 is the target admixed population

data <- "../qpgraph/plink_files/souilmi2024_bergstrom_plassias_shyam2025_transversions_for_qpadm_qpgraph"


test_for_hwd <- qpdstat(data, pop1="CoyoteCalifornia", pop2 = "HighlandWildDog", pop3 = "Dingo_Ancient_Nullarbor", pop4=masked_wgs_samples)
test_for_papua <- qpdstat(data, pop1="CoyoteCalifornia", pop2 = "Papua_Indonesia_518bp", pop3 = "Dingo_Ancient_Nullarbor", pop4=masked_wgs_samples)
test_for_ngsd <- qpdstat(data, pop1="CoyoteCalifornia", pop2 = "NewGuineaSingingDog", pop3 = "Dingo_Ancient_Nullarbor", pop4=masked_wgs_samples)

masked_wgs_samples[str_detect(masked_wgs_samples, "Curracurrang")]
masked_wgs_samples[str_detect(masked_wgs_samples, "bp")]

qpwave_ngsd_nullarbor <- qpwave(data, c("HighlandWildDog", "Papua_Indonesia_518bp", "NewGuineaSingingDog"), c("Dingo_Ancient_Nullarbor", masked_wgs_samples[str_detect(masked_wgs_samples, "Nullarbor")]),
                                verbose = TRUE)
qpwave_ngsd_se <- qpwave(data, c(masked_wgs_samples[str_detect(masked_wgs_samples, "Curracurrang")]), c("CoyoteCalifornia", "Dingo_Ancient_Nullarbor", "Papua_Indonesia_518bp", "Russia_Baikal_7k", ""),
                                verbose = TRUE)

source_pops <- c("Dingo_Ancient_Foreshore_697bp",
          "Dingo_Ancient_VIC_06", 
          "Dingo_Ancient_VIC_08", 
          "Dingo_Ancient_Curracurrang_D10",
          "Dingo_Ancient_Curracurrang_D13",
          "Dingo_Ancient_Curracurrang_D16",
          "Dingo_Ancient_Foreshore_715bp", 
          "Dingo_Ancient_Gambier_724bp", 
          "Dingo_Ancient_Gatton_13", 
          "Dingo_Ancient_NSW_SE")

modern_pops <- c("Dingo_FraserIsland_D02",
                 "Dingo_FraserIsland_E56",
                 "Dingo_Northern_Queensland_D04",
                 "Dingo_Northern_Queensland_D10",
                 "Dingo_FraserIsland_Y07",
                 "Dingo_FraserIsland_Y47",
                 "Dingo_Alpine_722g",
                 "Dingo_Alpine_Cooinda",
                 "Dingo_Alpine_D01",
                 "Dingo_Alpine_D05",
                 "Dingo_Alpine_D07",
                 "Dingo_Desert_Sandy",
                 "Dingo_Gibson_Desert_D08",
                 "Dingo_Kimberley_D09",
                 "Dingo_Northern_Territory_D06",
                 "Dingo_Simpson_Desert_D03")

snp_array_pops <- fread("cairns_samples.list", header=F, col.names = c("id"))$id

test_pops <- c("NewGuineaSingingDog", "Papua_Indonesia_518bp")

test_pops <- fread("test_pops.list", header = F, col.names = c("pop"))$pop

#curr_results <- qpf4ratio(data, c("CoyoteCalifornia","Dingo_Ancient_Nullarbor","Dingo_Ancient_VIC_08","Papua_Indonesia_518bp",'Dingo_Ancient_Curracurrang'))

qpf4ratio(data, c("AndeanFox","Zhokhov9500BP.CGG6","Dingo_FraserIsland_D02", "Dingo_Ancient_Nullarbor", "Papua_Indonesia_518bp"), verbose=T)

#filter_combined_results_archaic_deni <- filter(combined_results_archaic_deni, alpha < 1, alpha > 0, z > 2) 
combined_results <- data.frame()

for (source in snp_array_pops) {
  curr_results <- qpf4ratio(data, c('CoyoteCalifornia','NewGuineaSingingDog',source, "Dingo_Ancient_Nullarbor",  "Papua_Indonesia_518bp"), verbose=T)
  #curr_results <- qpf4ratio(data, c('CoyoteCalifornia','Dingo_Ancient_Nullarbor','Papua_Indonesia_518bp', "NewGuineaSingingDog",  source), verbose=T)
  print(curr_results)
  combined_results <- rbind(combined_results, curr_results)
}

#combined_results_filtered <- filter(combined_results, alpha < 1, alpha > 0, z > 2) 

write.table(combined_results,file="f4_ratio_Coyote_nulla_papua_ngsd_cairnsmodernx.csv")

cnulpapkarx <- fread("f4_ratio_Coyote_nulla_papua_karelia_easternsamps.csv")
cnulpapngsdx <- fread("f4_ratio_Coyote_nulla_papua_ngsd_easternsamps.csv")
cngsdxnullpap <- fread("f4_ratio_Coyote_ngsd_x_nulla_papua.csv")

##################
#### Karelia Dingo into Papua

cnulpapkarx %>% 
  left_join(., mds_n_metadata_extra_pop, by=join_by("pop5"=="rn")) %>%
  arrange(alpha) %>%
  mutate(pop5 = factor(pop5, levels=pop5)) %>%
  ggplot(aes(x=alpha, y=pop5)) +
  geom_errorbarh(aes(xmin=alpha-se, xmax=alpha+se), height=.2) +
  geom_point(
    aes(shape=type, fill=population), size=3
    ) +
  theme_classic()


##### dingo into Papua with ngsd


dingo_ancestry_in_ancient_ngsd <- cnulpapngsdx %>%
  filter(!str_detect(pop5, "Foreshore")) %>% 
  left_join(., mds_n_metadata_extra_pop, by=join_by("pop5"=="rn")) %>%
  arrange(alpha) %>%
  mutate(pop5 = factor(pop5, levels=pop5)) %>%
  ggplot(aes(x=alpha, y=pop5)) +
  geom_errorbarh(aes(xmin=alpha-se, xmax=alpha+se, size=as.factor(z > 2)), height=.1) +
  geom_point(
    aes(shape=type, fill=population), size=4
  ) +
  scale_fill_manual(name="Population", values=dingo_ngsd_pops_colours) +
  scale_size_manual(name="|Z| > 2", values=c("TRUE"=0.5, "FALSE"=0.1)) + 
  scale_shape_manual(name="Type", values=ancient_modern_shape) +
  labs(x="alpha", y="Source Sample") + 
  guides(fill = guide_legend(override.aes = list(shape = 21, size=3)), shape="none") + 
  theme_classic()
dingo_ancestry_in_ancient_ngsd
#### Papua into dingo with ngsd
anc_ngsd_ancestry <- cngsdxnullpap %>% 
  filter(!str_detect(pop3, "Foreshore")) %>% 
  left_join(., mds_n_metadata_extra_pop, by=join_by("pop3"=="rn")) %>%
  arrange(alpha) %>%
  mutate(pop3 = factor(pop3, levels=pop3)) %>%
  ggplot(aes(x=alpha, y=pop3)) +
  geom_errorbarh(aes(xmin=alpha-se, xmax=alpha+se, size=as.factor(z > 2)), height=.1) +
  geom_point(
    aes(shape=type, fill=population), size=4
  ) +
  scale_fill_manual(name="Population", values=dingo_ngsd_pops_colours) +
  scale_size_manual(name="|Z| > 2", values=c("TRUE"=0.5, "FALSE"=0.1)) + 
  scale_shape_manual(name="Type", values=ancient_modern_shape) +
  labs(x="alpha", y="Target Sample") + 
  guides(fill = "none", shape="none", size="none") + 
  theme_classic()
  
dingo_ancestry_in_ancient_ngsd + anc_ngsd_ancestry


dingo_into_ngsd <- cnulpapngsdx %>%
  filter(!str_detect(pop5, "Foreshore")) %>% 
  dplyr::select(pop5, alpha, se, z) %>%
  dplyr::rename(sample = pop5) %>%
  left_join(., mds_n_metadata_extra_pop, by=join_by("sample"=="rn")) %>%
  left_join(., metadata, by=join_by("sample"=="plink_name")) %>%
  left_join(., metadata_old_samples, by=join_by("sample"=="plink_name")) %>%
  mutate(sample_id_in_analysis := case_when(
    str_detect(sample, "Curracurrang") ~ paste0("Curracurrang.", date.y, "bp"),
    str_detect(sample, "Nullarbor") ~ paste0("Nullarbor.", date.y, "bp"),
    str_detect(sample, "FraserIsland") ~ str_replace(sample, "Dingo_FraserIsland", "Kgari"),
    is.na(sample_id_in_analysis) & str_detect(sample, "Dingo_") ~ str_replace(sample, "Dingo_", ""),
    str_detect(sample, "NewGuinea") ~ sample,
    .default = sample_id_in_analysis
  ),
  date := case_when(
    !is.na(date.x) ~ date.x,
    !is.na(date.y) ~ date.y,
    sample_id_in_analysis == "Gatton.nd" ~ 2500,
    .default = 0
  )) %>%
  mutate(test = "Dingo ancestry in ancient NGSD")

ngsd_into_dingo <- cngsdxnullpap %>%
  filter(!str_detect(pop3, "Foreshore")) %>% 
  dplyr::select(pop3, alpha, se, z) %>%
  dplyr::rename(sample = pop3) %>%
  left_join(., mds_n_metadata_extra_pop, by=join_by("sample"=="rn")) %>%
  left_join(., metadata, by=join_by("sample"=="plink_name")) %>%
  left_join(., metadata_old_samples, by=join_by("sample"=="plink_name")) %>%
  mutate(sample_id_in_analysis := case_when(
    str_detect(sample, "Curracurrang") ~ paste0("Curracurrang.", date.y, "bp"),
    str_detect(sample, "Nullarbor") ~ paste0("Nullarbor.", date.y, "bp"),
    str_detect(sample, "FraserIsland") ~ str_replace(sample, "Dingo_FraserIsland", "Kgari"),
    is.na(sample_id_in_analysis) & str_detect(sample, "Dingo_") ~ str_replace(sample, "Dingo_", ""),
    str_detect(sample, "NewGuinea") ~ sample,
    .default = sample_id_in_analysis
  ),
  date := case_when(
    !is.na(date.x) ~ date.x,
    !is.na(date.y) ~ date.y,
    sample_id_in_analysis == "Gatton.nd" ~ 2500,
    .default = 0
  )) %>%
  mutate(test = "NGSD ancestry in ancient dingoes")


ngsd_dingo_gene_flow <- rbind(dingo_into_ngsd, ngsd_into_dingo) %>%
  arrange(population, alpha) %>%
  mutate(sample_id_in_analysis = factor(sample_id_in_analysis, levels=unique(sample_id_in_analysis))) %>%
  ggplot(aes(x=alpha, y=sample_id_in_analysis)) +
  geom_errorbarh(aes(xmin=alpha-se, xmax=alpha+se, size=as.factor(z > 2)), height=.1) +
  geom_point(
    aes(fill=population), size=4, shape=21
  ) +
  facet_grid(~test) + 
  scale_fill_manual(name="Population", values=dingo_ngsd_pops_colours) +
  scale_size_manual(name="|Z| > 2", values=c("TRUE"=0.5, "FALSE"=0.1)) + 
  scale_shape_manual(name="Type", values=ancient_modern_shape) +
  labs(x="alpha", y="X") + 
  #guides(fill = "none", shape="none", size="none") + 
  theme_bw() +
  theme(
    panel.grid = element_blank()
  )


ggsave("ngsd_dingo_gene_flow.png", ngsd_dingo_gene_flow, dpi=500, height=4, width=8)

  
dingo_into_ngsd %>%
  filter(z > 2) %>%
  mutate(mean = mean(alpha)) %>%
  group_by(population) %>%
  mutate(pop_mean=mean(alpha)) %>%
  dplyr::select(mean, pop_mean)

ngsd_into_dingo %>%
  filter(z > 2) %>%
  mutate(mean = mean(alpha)) %>%
  group_by(population) %>%
  mutate(pop_mean=mean(alpha)) %>%
  dplyr::select(mean, pop_mean) %>%
  unique()



########## Masked f4 ratio
data <- "merge_masked_all_v2"
# Do Coyote and Nullarbor form a clade relative to HWD and papua?
qpdstat(data, pop1="AndeanFox", pop2 = "Dingo_Ancient_Nullarbor", pop3 = "HighlandWildDog", pop4="Papua_Indonesia_518bp",
        allsnps=TRUE, auto_only=F, f4mode=F)
qpdstat(data, pop1="AndeanFox", pop2 = "Dingo_Ancient_Nullarbor", pop3 = "HighlandWildDog", pop4="Papua_Indonesia_518bp",
        allsnps=TRUE, auto_only=F, f4mode=T)
qpdstat(data, pop1="CoyoteCalifornia", pop2 = "Papua_Indonesia_518bp", pop3 = "HighlandWildDog", pop4="Dingo_Ancient_Nullarbor",
        allsnps=TRUE, auto_only=F, f4mode=T)
qpdstat(data, pop1="CoyoteCalifornia", pop2 = "Papua_Indonesia_518bp", pop3 = "HighlandWildDog", pop4="Dingo_Ancient_Nullarbor",
        allsnps=TRUE, auto_only=F, f4mode=F)
qpdstat(data, pop1="CoyoteCalifornia", pop2 = "Dingo_Ancient_Nullarbor", pop4 = "NewGuineaSingingDog", pop3="Papua_Indonesia_518bp",
        allsnps=TRUE, auto_only=F, f4mode=F)
qpdstat(data, pop1="CoyoteCalifornia", pop2 = "Dingo_Ancient_Nullarbor", pop4 = "NewGuineaSingingDog", pop3="Papua_Indonesia_518bp",
        allsnps=TRUE, auto_only=F, f4mode=T)
masked_wgs_samples <- fread("qpadm_samples.list", header=F, col.names = c("id"))$id
combined_results <- data.frame()
source <- "Dingo_Ancient_Curracurrang_D16"



for (source in masked_wgs_samples) {
  curr_results <- qpf4ratio(curr_f4, c('CoyoteCalifornia','HighlandWildDog',source, "Papua_Indonesia_518bp", "Dingo_Ancient_Nullarbor"), verbose=T)
  #curr_results <- qpf4ratio(data, c('CoyoteCalifornia','Dingo_Ancient_Nullarbor','Papua_Indonesia_518bp', "HighlandWildDog",  source), verbose=T)
  print(curr_results)
  combined_results <- rbind(combined_results, curr_results)
}

#combined_results_filtered <- filter(combined_results, alpha < 1, alpha > 0, z > 2) 

write.table(combined_results,file="f4_ratio_Coyote_hwd_maskedmodernx_nulla_papua_2.csv")
#write.table(combined_results,file="f4_ratio_Coyote_nulla_papua_hwd_maskedmodernx.csv")

#f4_ratio_Coyote_hwd_maskedmodernx_nulla_papua <- combined_results

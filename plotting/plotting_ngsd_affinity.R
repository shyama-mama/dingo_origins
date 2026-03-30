pacman::p_load(
  dplyr,
  ggplot2,
  tidyr,
  stringr,
  rnaturalearth,
  rnaturalearthdata,
  sf
)


f4_matrix_treetest <- fread("~/Library/CloudStorage/Box-Box/projects/dingo/new_dingo/outgroup_f3/f4_matrix_treetest.csv")


all_sample_coords <- fread("~/Library/CloudStorage/Box-Box/projects/dingo/new_dingo/data/all_sample_coords.tsv")
all_sample_metadata <- fread("~/Library/CloudStorage/Box-Box/projects/dingo/new_dingo/metadata/metadata_for_wgs_samples.tsv")

f4_ratio_hwd_v3 <- fread("~/Library/CloudStorage/Box-Box/projects/dingo/new_dingo/data/f4_ratio_Coyote_hwd_x_nulla_papua_v3.csv") %>%
  dplyr::select(pop3, alpha, se, z)
f4_ratio_hwd_v2 <- fread("~/Library/CloudStorage/Box-Box/projects/dingo/new_dingo/data/f4_ratio_Coyote_hwd_x_nulla_papua.csv") %>%
  dplyr::select(pop3, alpha, se, z)
f4_ratio_ngsd_v3 <- fread("~/Library/CloudStorage/Box-Box/projects/dingo/new_dingo/data/f4_ratio_Coyote_ngsd_x_nulla_papua_v3.csv") %>%
  dplyr::select(pop3, alpha, se, z)
merge(f4_ratio_hwd_v2, f4_ratio_hwd_v3, by="pop3") %>%
  filter(alpha.x != alpha.y) %>%
  ggplot(aes(x=alpha.x, y=alpha.y)) +
  geom_abline(intercept = 0) + 
  geom_point() +
  geom_errorbar(aes(ymin=alpha.y-se.y, ymax=alpha.y+se.y)) +
  geom_errorbarh(aes(xmin=alpha.x-se.x, xmax=alpha.x+se.x))
  

f4_ratio_all_hwd <- left_join(f4_ratio_hwd_v3, all_sample_coords, by=join_by("pop3"=="sample"))%>%
  left_join(., all_sample_metadata, by=join_by("pop3"=="sample"))

f4_ratio_all_hwd_ngsd <- left_join(f4_ratio_all_hwd, f4_ratio_ngsd_v3, by="pop3")

f4_ratio_all_hwd_ngsd %>%
  ggplot(aes(x=alpha.x, y=alpha.y)) +
  geom_abline(intercept = 0) + 
  geom_point() +
  geom_errorbar(aes(ymin=alpha.y-se.y, ymax=alpha.y+se.y)) +
  geom_errorbarh(aes(xmin=alpha.x-se.x, xmax=alpha.x+se.x))

ggplotly(f4_ratio_all_hwd %>%
  ggplot(aes(x=longitude, y=alpha, label=pop3, colour=population_label)) +
  geom_point())

ggplotly(f4_ratio_all_hwd %>%
           filter(population_label != "") %>%
           ggplot(aes(x=population_label, y=alpha, label=pop3, colour=population_label)) +
           geom_violin(scale = "width")+
           geom_point())


qpadm_pap_df <- qpadm_hwd_pap_ranked_best %>%
  dplyr::select(target, Dingo_Ancient_Nullarbor, Papua_Indonesia_518bp, p, se_Dingo_Ancient_Nullarbor, se_Papua_Indonesia_518bp, z_Dingo_Ancient_Nullarbor, z_HighlandWildDog, z_Papua_Indonesia_518bp) %>%
  replace(is.na(.), 0) 

f4_ratio_all_hwd_qpadm <- left_join(f4_ratio_all_hwd, qpadm_pap_df, by=join_by("pop3"=="target"))



australia <- ne_countries(scale = "large", country = "Australia", returnclass = "sf")

qpadm_ngsd_plot <- ggplot() +
  geom_sf(data=australia, fill=NA)+
  geom_point(data=f4_ratio_all_hwd_qpadm %>% filter(!is.na(Papua_Indonesia_518bp), age=="Modern", population_label != "Captive"), aes(x=longitude, y=latitude, fill=Papua_Indonesia_518bp), shape=21, size=3) +
  #geom_point(data=f4_ratio_all_hwd_qpadm %>% filter(!is.na(Papua_Indonesia_518bp), age=="Ancient"), aes(x=longitude, y=latitude, fill=Papua_Indonesia_518bp, label=pop3), shape=24, size=4) +
  scale_fill_continuous("\nNew Guinean Ancestry", low="white", high="black", limits=c(0, 0.26), breaks=c(0.05, 0.15, 0.25)) +
  #scale_colour_manual(values=c("white", "black")) + 
  theme_void() +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, direction = "horizontal", text.size=5)) +
  coord_sf(ylim=c(-42,-7))+
  theme(
    legend.position = c(0.3,0.1),
    legend.text = element_text(angle = 45, size=9)
  )





qpadm_per_population <- f4_ratio_all_hwd_qpadm %>%
  mutate(
    Papua_Indonesia_518bp := case_when(
      pop3 == "Dingo_Ancient_Curracurrang_D16" ~ 0.4085631,
      .default = Papua_Indonesia_518bp
    ),
    se_Papua_Indonesia_518bp := case_when(
      pop3 == "Dingo_Ancient_Curracurrang_D16" ~ 0.03145031,
      .default = se_Papua_Indonesia_518bp
    )
  ) %>%
  filter(!is.na(Papua_Indonesia_518bp), !population_label %in% c("Captive", "", "North")) %>%
  mutate(population_label = factor(population_label, levels=c("Ancient West", "West", "North", "Central", "Mallee", "Alpine", "East", "K'gari", "HWD", "NGSD", "Ancient NG", "Ancient South", "Ancient East")),) %>%
ggplot(aes(label=pop3)) +
  geom_violin(aes(x=population_label, y=Papua_Indonesia_518bp, fill=population_label), colour=NA, alpha=0.4, scale="width") + 
  geom_pointrange(aes(x=population_label, y=Papua_Indonesia_518bp, ymin = Papua_Indonesia_518bp-se_Papua_Indonesia_518bp, ymax = Papua_Indonesia_518bp+se_Papua_Indonesia_518bp, shape=age, fill=population_label), 
                  position=position_jitter(width=0.1), 
                  linetype='dotted', stroke=0.5, size=0.5, alpha=0.7) +
  #geom_point(aes(x=population_label, y=Papua_Indonesia_518bp, shape=age, fill=population_label), size=3) +
  scale_shape_manual(values=c("Modern"=21, "Ancient"=24)) + 
  labs(x="Population", y="New Guinean Ancestry") + 
  guides(shape="none", fill="none", colour="none") +
  scale_fill_manual(values=pop_colours) + 
  #scale_colour_manual(values=pop_colours) + 
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle=90),
  )

ggplotly(qpadm_per_population)

figure_qpadm_ngsd_ancestry <- qpadm_per_population | qpadm_ngsd_plot
ggsave("figure_qpadm_ngsd_ancestry.png", figure_qpadm_ngsd_ancestry, dpi=300, height=5, width=8)


ggplotly(f4_ratio_all_hwd_qpadm %>% 
  filter(!is.na(Papua_Indonesia_518bp), age=="Modern") %>%
  ggplot(aes(x=longitude, y=Papua_Indonesia_518bp, label=pop3, colour=population_label)) +
  geom_point() +
  theme_bw())

ggplotly(f4_ratio_all_hwd_qpadm %>% 
           filter(!is.na(Papua_Indonesia_518bp), age=="Modern") %>%
           ggplot(aes(x=latitude, y=Papua_Indonesia_518bp, colour=population_label)) +
           geom_point() +
           theme_bw())



ggplot() +
  geom_sf(data=australia)+
  geom_point(data=f4_ratio_all_hwd %>% mutate(alpha := case_when(alpha < 0 ~ 0, .default = alpha)) %>% filter(age == "Modern", z < 3), aes(x=longitude, y=latitude, fill=alpha, colour=as.factor(z>3)), shape=21, size=2) +
  geom_point(data=f4_ratio_all_hwd %>% mutate(alpha := case_when(alpha < 0 ~ 0, .default = alpha)) %>% filter(age == "Modern", z >= 3), aes(x=longitude, y=latitude, fill=alpha, colour=as.factor(z>3)), shape=21, size=3) +
  scale_fill_continuous(low="white", high="black") +
  scale_colour_manual(values=c("white", "black")) + 
  theme_void() +
  coord_sf(ylim=c(-42,-5))


ggplotly(ggplot() +
  geom_sf(data=australia)+
  geom_point(data=f4_ratio_all_hwd, aes(x=longitude, y=latitude, fill=population_label,label=pop3), shape=21, size=2) +
  theme_void() +
  coord_sf(ylim=c(-42,-5)))

left_join(f4_ratio_all_hwd, f4_summary_for_all, by=join_by("pop3"=="sample")) %>%
  ggplot() +
  geom_point(aes(x=affinity, y=alpha))



ggplot() +
  geom_sf(data=australia)+
  geom_point(data=f4_summary_for_all %>% filter(age == "Modern"), aes(x=longitude, y=latitude, fill=affinity), shape=21, size=3) +
  scale_fill_continuous(low="white", high="black") +
  scale_colour_manual(values=c("white", "black")) + 
  theme_void() +
  coord_sf(ylim=c(-42,-5))

merge_alpha_f4 <- merge(f4_ratio, modern_samples_f4_w_coords, by.x="pop3", by.y="pop4")
merge_alpha_f4 %>%
  ggplot(aes(x=alpha, y=est))+
  geom_smooth(method="lm") + 
  geom_errorbar(aes(ymin=est-se.y, ymax=est+se.y, colour=as.factor(z.y>2))) +
  geom_errorbarh(aes(xmin=alpha-se.x, xmax=alpha+se.x, colour=as.factor(z.x>2))) +
  geom_point() +
  theme_minimal()


merge_alpha_f4 %>%
  ggplot(aes(x=longitude, y=alpha))+
  geom_smooth(method="lm") + 
  geom_errorbar(aes(ymin=alpha-se.x, ymax=alpha+se.x, colour=as.factor(z.y>2))) +
  geom_point() +
  theme_minimal()


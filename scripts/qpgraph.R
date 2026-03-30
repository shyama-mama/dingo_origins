pacman::p_load(
  ggplot2,
  plotly,
  admixtools,
  dplyr,
  data.table,
  patchwork
)


setwd("/Users/shyamamama/Library/CloudStorage/Box-Box/projects/dingo/new_dingo/qpgraph/")


get_f2 <- function(eigen_prefix, pops, output_prefix) {
  extracted_f2 <- extract_f2(eigen_prefix, pops=pops, outdir=output_prefix,
                             adjust_pseudohaploid = TRUE, auto_only = FALSE, 
                             verbose = TRUE, overwrite = TRUE)
  
  f2 <- f2_from_precomp(output_prefix)
  
  return(f2)
}

top_find_graph <- function(f2, output_prefix, numadmix=0, outpop="Sheep") {
  graph <- find_graphs(
    f2, numadmix=numadmix,
    outpop = outpop,
    verbose = TRUE
  )
  
  subset <- graph %>% slice_min(score, n = 6)
  num_graphs <- nrow(subset)
  
  one <- subset$edges[[1]] %>% plot_graph(title = paste(round(subset$score[[1]], 2)))
  two <- subset$edges[[2]] %>% plot_graph(title = paste(round(subset$score[[2]], 2)))
  three <- subset$edges[[3]] %>% plot_graph(title = paste(round(subset$score[[3]], 2)))
  
  if (num_graphs < 6) {
    print("Getting top 3!")
    layout <- "
    AB
    C#
    "
    multiple_plots <- one + two + three + plot_layout(design = layout)
    
  } else {
    print("Getting top 6!")
    layout <- "
    ABC
    DEF
    "
    four <- subset$edges[[4]] %>% plot_graph(title = paste(round(subset$score[[4]], 2)))
    five <- subset$edges[[5]] %>% plot_graph(title = paste(round(subset$score[[5]], 2)))
    six <- subset$edges[[6]] %>% plot_graph(title = paste(round(subset$score[[6]], 2)))
    multiple_plots <- one + two + three + four + five + six + plot_layout(design = layout)
  }
  
  output_file <- paste0(output_prefix, "_", numadmix, "admix")
  ggsave(paste0(output_file, ".png"), multiple_plots, dpi=600, width=10, height=6)
  ggsave(paste0(output_file, ".pdf"), multiple_plots, dpi=600, width=10, height=6)
  
  
  return(list(plot=multiple_plots, top_graphs=subset))
}

get_residual <- function(f2, edge, f3basepop="Sheep", output_prefix="residual", numstart=100){
  qp_graph = qpgraph(f2, edge, return_pvalue = TRUE, return_fstats = TRUE, 
                     verbose = TRUE, allsnps = TRUE, boot = TRUE, numstart = numstart, f3basepop=f3basepop)
  
  annotation <- paste0("worst_residual: ", round(qp_graph$worst_residual, 2), "\nscore: ", round(qp_graph$score, 2)) 
  plotted_graph = qp_graph$edges %>% plot_graph(highlight_unidentifiable = FALSE, title = annotation)
  
  qp_f3_residuals <- as.data.table(qp_graph$f3)
  qp_f3_residuals[, Sig:=""]
  qp_f3_residuals[abs(z)>=2, Sig:="*"]
  qp_f3_residuals[abs(z)>=3, Sig:="**"]
  
  qp_f3_residuals[, pop2:=factor(pop2, levels=qp_f3_residuals[, .N, by=pop2][order(N), pop2])]
  qp_f3_residuals[, pop3:=factor(pop3, levels=qp_f3_residuals[, .N, by=pop3][order(N), pop3])]
  residual_plot <- ggplot(qp_f3_residuals, aes(x=pop2, y=pop3, fill=z)) +
    geom_tile(color="black", size=0.25) +
    geom_text(aes(label=round(z,2), size=2), show.legend=F) +
    scale_fill_gradient2(high = "red",
                         mid = "white ",
                         low = "blue") +
    labs(fill="Residual (Z)") +
    theme_bw() +
    theme(axis.title=element_blank(),
          axis.text.x=element_text(hjust=1, angle=45),
          legend.position=c(0.125, 0.25),
          legend.background=element_rect(color="black", size=0.5))
  
  layout <- "
  A
  B"
  
  merged_plot <- plotted_graph + residual_plot + plot_layout(design = layout)
  
  ggsave(paste0(output_prefix, ".png"), merged_plot, dpi=600, width=10, height=10)
  ggsave(paste0(output_prefix, ".pdf"), merged_plot, dpi=600, width=10, height=10)
  
  return(list(qp_graph=qp_graph, plot=merged_plot))
}

get_all_residual <- function(f2, graph_output, f3basepop="Sheep", output_prefix="residual") {
  num_graphs <- nrow(graph_output$top_graphs)
  if(num_graphs < 6) {
    total <- 3
  } else {
    total <- 6
  }
  
  for (i in seq(1,total)) {
    get_residual(f2, graph_output$top_graphs$edges[[i]], f3basepop=f3basepop, output_prefix=paste0(output_prefix, "_", i))
  }
}


## Find graphs
call_set <- "plink_files/merge_masked_all"
output_prefix <- "plink_files/Nullarbor_East_Papua_Hwd"
f2 <- get_f2(call_set, c("AndeanFox", "WolfSyria",  "Dingo_Ancient_East", "Dingo_Ancient_Nullarbor", "Israel_7000BP.THRZ02", "Karelia_Veretye", "Papua_Indonesia_518bp", "HighlandWildDog"),
             output_prefix)

# f2_only_old <- get_f2(call_set, c(), "plink_files/souilmi2024_bergstrom_plassias_shyam2025_transversions_for_qpadm_qpgraph_3")
# f2_old_new <- get_f2(call_set, c(), "plink_files/souilmi2024_bergstrom_plassias_shyam2025_transversions_for_qpadm_qpgraph_2")
# f2_only_new <- get_f2(call_set, c(), "plink_files/souilmi2024_bergstrom_plassias_shyam2025_transversions_for_qpadm_qpgraph_1")

graph_output_o <- top_find_graph(f2, "souilmi2024_east_hwd_new_papua_find_graph", numadmix=0, outpop="AndeanFox")
graph_output_o_1 <- top_find_graph(f2, "souilmi2024_east_hwd_new_papua_find_graph", numadmix=1, outpop="AndeanFox")
graph_output_o_2 <- top_find_graph(f2, "souilmi2024_east_hwd_new_papua_find_graph", numadmix=2, outpop="AndeanFox")
graph_output_o_3 <- top_find_graph(f2, "souilmi2024_east_hwd_new_papua_find_graph", numadmix=3, outpop="AndeanFox")
graph_output_o_4 <- top_find_graph(f2, "souilmi2024_east_hwd_new_papua_find_graph", numadmix=4, outpop="AndeanFox")


get_all_residual(f2, graph_output_o, f3basepop="AndeanFox", output_prefix="souilmi2024_plus_new_papua_find_graph_0")
get_all_residual(f2, graph_output_o_1, f3basepop="AndeanFox", output_prefix="souilmi2024_plus_new_papua_find_graph_1")
get_all_residual(f2, graph_output_o_2, f3basepop="AndeanFox", output_prefix="souilmi2024_plus_new_papua_find_graph_2")
get_all_residual(f2, graph_output_o_3, f3basepop="AndeanFox", output_prefix="souilmi2024_plus_new_papua_find_graph_3")
get_all_residual(f2, graph_output_o_4, f3basepop="AndeanFox", output_prefix="souilmi2024_plus_new_papua_find_graph_4")

get_all_residual(f2_only_old, graph_output_o, f3basepop="AndeanFox", output_prefix="souilmi2024_only_new_papua_find_graph_0")
get_all_residual(f2_only_old, graph_output_o_1, f3basepop="AndeanFox", output_prefix="souilmi2024_only_new_papua_find_graph_1")
get_all_residual(f2_only_old, graph_output_o_2, f3basepop="AndeanFox", output_prefix="souilmi2024_only_new_papua_find_graph_2")
get_all_residual(f2_only_old, graph_output_o_3, f3basepop="AndeanFox", output_prefix="souilmi2024_only_new_papua_find_graph_3")
get_all_residual(f2_only_old, graph_output_o_4, f3basepop="AndeanFox", output_prefix="souilmi2024_only_new_papua_find_graph_4")


### manual graph
call_set <- "plink_files/merge_masked_all"
manual_graph_1 <- fread("souilmi_manual_graph_1.tsv")
output_prefix <- "plink_files/Nullarbor_Curra_NGSD"
f2 <- get_f2(call_set, c("AndeanFox", "WolfSyria",  "Dingo_Ancient_Curracurrang", "Dingo_Ancient_Nullarbor", "Israel_7000BP.THRZ02", "Karelia_Veretye", "NewGuineaSingingDog"),
             output_prefix)

manual_graph_1 %>% plotly_graph()

prefix <- "souilmi2024_original_6M"
residual <- get_residual(f2, manual_graph_1, f3basepop="AndeanFox", output_prefix=paste0(prefix, "_manual_graph"), numstart=1000)


manual_graph_1 <- fread("souilmi_manual_graph_1_w_pap.tsv")
output_prefix <- "plink_files/Nullarbor_Curra_Pap"
f2 <- get_f2(call_set, c("AndeanFox", "WolfSyria",  "Dingo_Ancient_Curracurrang", "Dingo_Ancient_Nullarbor", "Israel_7000BP.THRZ02", "Karelia_Veretye", "Papua_Indonesia_518bp"),
             output_prefix)

manual_graph_1 %>% plotly_graph()

prefix <- "souilmi2024_6M_w_pap"
residual <- get_residual(f2, manual_graph_1, f3basepop="AndeanFox", output_prefix=paste0(prefix, "_manual_graph"), numstart=1000)

manual_graph_1 <- fread("souilmi_manual_graph_1_w_hwd.tsv")
output_prefix <- "plink_files/Nullarbor_Curra_Hwd"
f2 <- get_f2(call_set, c("AndeanFox", "WolfSyria",  "Dingo_Ancient_Curracurrang", "Dingo_Ancient_Nullarbor", "Israel_7000BP.THRZ02", "Karelia_Veretye", "HighlandWildDog"),
             output_prefix)

manual_graph_1 %>% plotly_graph()

prefix <- "souilmi2024_6M_w_hwd"
residual <- get_residual(f2, manual_graph_1, f3basepop="AndeanFox", output_prefix=paste0(prefix, "_manual_graph"), numstart=1000)

manual_graph_1 <- fread("souilmi_manual_graph_1_w_pap_w_east.tsv")
output_prefix <- "plink_files/Nullarbor_East_Pap"
f2 <- get_f2(call_set, c("AndeanFox", "WolfSyria",  "Dingo_Ancient_East", "Dingo_Ancient_Nullarbor", "Israel_7000BP.THRZ02", "Karelia_Veretye", "Papua_Indonesia_518bp"),
             output_prefix)

manual_graph_1 %>% plotly_graph()

prefix <- "souilmi2024_6M_w_pap_east"
residual <- get_residual(f2, manual_graph_1, f3basepop="AndeanFox", output_prefix=paste0(prefix, "_manual_graph"), numstart=1000)

manual_graph_1 <- fread("souilmi_manual_graph_1_w_pap_w_seast.tsv")
output_prefix <- "plink_files/Nullarbor_SEast_Pap"
call_set <- "plink_files/merge_masked_all_se"
f2 <- get_f2(call_set, c("AndeanFox", "WolfSyria",  "Dingo_Ancient_SEast", "Dingo_Ancient_Nullarbor", "Israel_7000BP.THRZ02", "Karelia_Veretye", "Papua_Indonesia_518bp"),
             output_prefix)

manual_graph_1 %>% plotly_graph()

prefix <- "souilmi2024_6M_w_pap_seast"
residual <- get_residual(f2, manual_graph_1, f3basepop="AndeanFox", output_prefix=paste0(prefix, "_manual_graph"), numstart=1000)

#### man graph 
manual_graph_2 <- fread("souilmi_manual_graph_1_w_pap.tsv")

manual_graph_2 %>% plotly_graph()

prefix <- "souilmi2024_only_new_papua_2"
residual <- get_residual(f2_only_old, manual_graph_2, f3basepop="AndeanFox", output_prefix=paste0(prefix, "_manual_graph"), numstart=1000)



############## 
manual_graph_1_org <- fread("souilmi_manual_graph_1.tsv")

manual_graph_1_org %>% plotly_graph()

prefix <- "souilmi2024_original_1"
residual <- get_residual(f2_only_new, manual_graph_1_org, f3basepop="AndeanFox", output_prefix=paste0(prefix, "_manual_graph"), numstart=1000)

manual_graph_2_org <- fread("souilmi_manual_graph_2.tsv")

manual_graph_2_org %>% plotly_graph()

prefix <- "souilmi2024_original_2"
residual <- get_residual(f2_only_new, manual_graph_2_org, f3basepop="AndeanFox", output_prefix=paste0(prefix, "_manual_graph"), numstart=1000)






########## 
manual_graph_1 <- fread("souilmi_manual_graph_1_w_pap_w_east.tsv")
output_prefix <- "plink_files/Nullarbor_East_Pap_19"
call_set <- "plink_files/souilmi2024_bergstrom_plassias_shyam2025_transversions_for_qpadm_qpgraph"
f2 <- get_f2(call_set, c("AndeanFox", "WolfSyria",  "Dingo_Ancient_East", "Dingo_Ancient_Nullarbor", "Israel_7000BP.THRZ02", "Karelia_Veretye", "Papua_Indonesia_518bp"),
             output_prefix)

manual_graph_1 %>% plotly_graph()

prefix <- "souilmi2024_19M_w_pap_w_east"
residual <- get_residual(f2, manual_graph_1, f3basepop="AndeanFox", output_prefix=paste0(prefix, "_manual_graph"), numstart=1000)

manual_graph_2 <- fread("souilmi_manual_graph_2_w_pap_w_east.tsv")
prefix <- "souilmi2024_19M_w_pap_w_east_2"
residual <- get_residual(f2, manual_graph_2, f3basepop="AndeanFox", output_prefix=paste0(prefix, "_manual_graph"), numstart=1000)

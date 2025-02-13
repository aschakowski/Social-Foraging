visualize_depletion_results <- function(){

comparison = readRDS(file = "utils/data/processed_data/model_fit/depletion_models/comparison/depletion_comparison.rds")
comparison = as.data.frame(comparison)
comparison$name = c("Exponential", "Holling III", "Michaelis-\nMenten", "Linear")
comparison$model_id = 1:4
comparison_depletion = comparison %>%
  ggplot(aes(x = elpd_diff, y = fct_reorder(name, model_id))) + 
  geom_col(color = "black", width = .6) +
  geom_errorbarh(aes(xmin = elpd_diff - 4*se_diff, xmax = elpd_diff + 4*se_diff), height = .3, size = .8) + 
  ggtitle("Model Comparison Patch Depletion") + 
  xlab(expression(Delta~"ELPD (rel. to best)")) + 
  ylab("Model") + 
  #labs(tag = "a")+
  geom_vline(xintercept = 0)+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7, angle = 65, hjust = .5),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.9, .2),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_rect(colour = "white", fill = "white"),
        text = element_text(family = "sans"))
comparison_depletion

ggsave(comparison_depletion, file = "4_incremental_strategies/output/depletion_comparison.png", 
       width = 8, height = 8, bg = "white", dpi = 1200, units = "cm")
ggsave(comparison_depletion, file = "4_incremental_strategies/output/depletion_comparison.svg", 
       width = 8, height = 8, bg = "white", dpi = 1200, units = "cm")
ggsave(comparison_depletion, file = "4_incremental_strategies/output/depletion_comparison.pdf", 
       width = 8, height = 8, bg = "white", dpi = 1200, units = "cm")
ggsave(comparison_depletion, file = "4_incremental_strategies/output/depletion_comparison.tiff", 
       width = 8, height = 8, bg = "white", dpi = 1200, units = "cm")

}
################################################################################
# END
################################################################################

library(tidyverse)
library(colorspace)
library(ggrepel)

dark_okabe <- darken(c("#E69F00", "#56B4E9", "grey60"), amount = 0.2) 

dir.create("./plots")

pca_me3 <- read.table("./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/bigwigs/me3_pca.tab", header = TRUE) %>%
  dplyr::mutate(variance = 100*Eigenvalue/sum(Eigenvalue)) %>%
  dplyr::filter(Component <= 2) %>%
  mutate(Component = paste0("PC", Component, " (", round(variance, digits = 1), "%)")) %>%
  dplyr::select(-c("Eigenvalue", "variance")) %>%
  set_names("Component", "WT me3 1", "WT me3 2", 
            "C5 me3 1", "C5 me3 2", 
            "C9 me3 1", "C9 me3 2") %>%
  pivot_longer(-Component, values_to = "value", names_to = "Sample") %>%
  pivot_wider(names_from = Component) %>%
  mutate(Condition = sub("(me3).*", "\\1", Sample)) %>%
  mutate(Condition = case_when(Condition == "C5 me3" ~ "C5 H3K27me3",
                               Condition == "C9 me3" ~ "C9 H3K27me3",
                               Condition == "WT me3" ~ "WT H3K27me3"))

me3_plot_pca <- pca_me3 %>% 
  ggplot(aes(x = `PC1 (59.9%)`, y = `PC2 (29.3%)`, colour = Condition, fill = Condition, label = Sample)) +
  geom_point(size=6, stroke = 1, shape = 21, alpha = 0.5) +
  geom_text_repel(nudge_x = 0.05, nudge_y = 0.15, point.padding = 1, show.legend = FALSE) +
  theme_classic(base_size = 22) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#000000")) + 
  scale_colour_manual(values = dark_okabe)

ggsave(file = "./plots/pca_me3.pdf", plot = me3_plot_pca, device = cairo_pdf(), 
       width = 8, height = 3.5, dpi = 800)

pca_ac <- read.table("./CUT_and_RUN/1_0_Snakemake_run_sep_replicates/bigwigs/ac_pca.tab", header = TRUE) %>%
  dplyr::mutate(variance = 100*Eigenvalue/sum(Eigenvalue)) %>%
  dplyr::filter(Component <= 2) %>%
  mutate(Component = paste0("PC", Component, " (", round(variance, digits = 1), "%)")) %>%
  dplyr::select(-c("Eigenvalue", "variance")) %>%
  set_names("Component", "WT ac 1", "WT ac 2", 
            "C5 ac 1", "C5 ac 2", 
            "C9 ac 1", "C9 ac 2") %>%
  pivot_longer(-Component, values_to = "value", names_to = "Sample") %>%
  pivot_wider(names_from = Component) %>%
  mutate(Condition = sub("(ac).*", "\\1", Sample)) %>%
  mutate(Condition = case_when(Condition == "C5 ac" ~ "C5 H3K27ac",
                               Condition == "C9 ac" ~ "C9 H3K27ac",
                               Condition == "WT ac" ~ "WT H3K27ac"))


ac_plot_pca <- pca_ac %>% 
  ggplot(aes(x = `PC1 (85.2%)`, y = `PC2 (11.8%)`, colour = Condition, fill = Condition, label = Sample)) +
  geom_point(size=6, stroke = 1, shape = 21, alpha = 0.5) +
  geom_text_repel(point.padding = 1, show.legend = FALSE) +
  #coord_fixed() +
  expand_limits(y = c(-0.65, 0.45)) +
  theme_classic(base_size = 22) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#000000")) + 
  scale_colour_manual(values = dark_okabe)


ggsave(file = "./plots/pca_ac.pdf", plot = ac_plot_pca, device = cairo_pdf(), 
       width = 10, height = 3, dpi = 800)
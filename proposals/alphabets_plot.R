library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)

size_mod <- 1

my_theme <- theme(
  axis.text = element_text(size=13 + size_mod),
  axis.title.x = element_text(size=14 + size_mod, vjust = -1),
  axis.title.y = element_text(size=14 + size_mod, vjust = 1),
  
  legend.background = element_rect(fill = "NA"),
  legend.key = element_rect(fill = "NA", color = "NA", size = 0.5),
  legend.position = "bottom",
  #uncomment for publications
  legend.key.size = unit(0.1, "inches"),
  #legend.margin = unit(-0.25, "lines"),
  legend.text = element_text(size=13 + size_mod), 
  legend.title = element_text(size=15 + size_mod),
  
  panel.grid.major = element_line(color="grey", linetype = "dashed", size = 0.5),
  panel.grid.major = element_line(color="lightgrey", 
                                  linetype = "dashed", size = 0.5),
  panel.background = element_rect(fill = "transparent", color = "black"),
  
  plot.background=element_rect(fill = "transparent",
                               color = "transparent"),
  #uncomment for publications
  plot.margin = unit(rep(0.2, 4), "inches"),
  plot.title = element_text(size=20 + size_mod),
  
  strip.background = element_rect(fill = "NA", color = "NA"),
  strip.text = element_text(size=13 + size_mod, face = "bold")
)

freq_nondeg <- read.csv2("./proposals/freq_nondeg.csv")[, -1]
freq_deg <- read.csv2("./proposals/freq_deg.csv")[, -1]

do_pca <- function(x) 
  x %>% 
  select(-type, -taxon) %>% 
  prcomp(center = TRUE, scale = TRUE) %>% 
  getElement("x") %>% 
  data.frame() %>% 
  select(1, 2) %>% 
  cbind(select(x, type, taxon), .) %>% 
  mutate(type_nice = factor(type, labels = c("dojrzałe białko\n", "peptyd sygnałowy\n")),
         taxon_nice = factor(taxon, labels = c("inne", "Plasmodium"))) %>% 
  mutate(both = paste0(type_nice, "(", taxon_nice, ")\n"))

dat_deg <- do_pca(freq_deg) 
dat_nondeg <- do_pca(freq_nondeg)

# dat_deg <- read.table("PCAgr.txt", header = TRUE, sep = "\t")
# dat_nondeg <- read.table("PCA.txt", header = TRUE, sep = "\t")
# colnames(dat_deg) <- c("both", "PC1", "PC2")
# colnames(dat_nondeg) <- c("both", "PC1", "PC2")

plot_pca <- function(x)
  ggplot(x, aes(x = PC1, y = PC2, color = both, shape = both, fill = both)) + 
  geom_density_2d(color = "black", contour = TRUE) +
  #geom_point() +
  stat_density2d(aes(fill=both,alpha=..level..), color = "black", contour = TRUE, geom="polygon") +
  scale_linetype_discrete("") +
  scale_fill_discrete("") +
  scale_shape_discrete("") +
  scale_color_discrete("") +
  scale_alpha_continuous(range = c(0.25, 0.4)) +
  guides(alpha = FALSE) +
  my_theme

p_nondeg <- plot_pca(dat_nondeg) + guides(fill = FALSE, linetype = FALSE) + ggtitle("Pełny alfabet")
p_deg <- plot_pca(dat_deg) + ggtitle("Skrócony alfabet signalHsmm")

cairo_ps("./proposals/signalHsmm.eps", width = 9, height = 9, onefile = FALSE)
grid.arrange(p_nondeg, p_deg, nrow = 2, ncol = 1, heights = c(0.5, 0.5))
dev.off()
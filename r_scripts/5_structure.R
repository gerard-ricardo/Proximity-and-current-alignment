# structure analysis (dart - working) ----------------------------------------------------
##notes
# data wrangling ----------------------------------------------------------
# structure ---------------------------------------------------------------
tic("Running structure analysis")
struct_adult = gl.run.structure(data_gl_adult_unique, verbose = 3, burnin = 2000, numreps = 4000, k.range = 1:5, num.k.rep = 2, 
                                seed = 1, noadmix=FALSE, exec = "C:/Users/gerar/OneDrive/Documents/structure/structure.exe")
toc() 
tic("Running structure analysis")
struct_adult = gl.run.structure(data_gl_adult_unique, verbose = 3, burnin = 8000, numreps = 20000, k.range = 2:4, num.k.rep = 2,
                                seed = 1, noadmix = FALSE, exec = "C:/Users/gerar/OneDrive/Documents/structure/structure.exe")
toc()
load("./Rdata/struct_adult_test.RData")
str(struct_adult)
ev <- gl.evanno(struct_adult, plot.out = TRUE)
gl.evanno(struct_adult)
data_gl_adult_unique1 = data_gl_adult_unique
pop(data_gl_adult_unique1) <- factor(rep("MergedPop", length(indNames(data_gl_adult_unique1))))
gl.report.hwe(data_gl_adult_unique1)
gl.filter.hwe(data_gl_adult_unique1, alpha = 0.05, mult.comp.adj = TRUE, mult.comp.adj.method = "fdr")
qmat <- dartR::gl.plot.structure(struct_adult, K = 3, colors_clusters = list("dodgerblue", "mediumseagreen", "salmon", 'pink'), clumpak = T, save2tmp = T)
qmat <- dartR::gl.plot.structure(struct_adult, K = 2, colors_clusters = list("dodgerblue", "mediumseagreen"), clumpak = T, save2tmp = F, 
                                 plot.out = F)
head(qmat)
p3 = gl.print.reports(1)
gl.list.reports()
Q_melt <- do.call("rbind", lapply(qmat, reshape2::melt, id.vars = c("Label", "K", "orig.pop", "ord"), variable.name = "Cluster" ))
Q_melt$orig.pop <-
  factor(Q_melt$orig.pop, levels = unique(struct_adult[[1]]$q.mat$orig.pop))
p3 <- ggplot(Q_melt, aes_(x= ~ factor(ord), y = ~value, fill = ~Cluster)) +
  geom_col(color = "black", size = 0.25, width = 1) +
  facet_grid(K ~ orig.pop , scales = "free", space = "free") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(
    breaks = unique(Q_melt$ord), labels = unique(Q_melt$Label), expand = c(0, 0)) +
  scale_fill_manual(values = c("dodgerblue", "mediumseagreen", "salmon")) +
  theme_sleek2() +
  theme(panel.spacing = unit(0, "lines"), panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12, angle = 90),
        axis.text.x = element_text(size = 8,angle = 90,vjust = 0.5,hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank() ,
        legend.position = "none") +
  labs(x = 'Individuals', y = 'K=3')
p3
load("./Rdata/structure_plot.RData")
gl.map.structure(qmat = qmat, x = data_gl_filtered_adult, K = 3, scalex = 1, scaley = 0.5)

class(pairwise_dist_plot)
class(polar_plot)
class(arranged_plots)
class(sankey_plot)
class(distance_map_plot)
library(gridExtra)
library(grid)
library(plotrix)
library(patchwork)
library(cowplot)
library(htmlwidgets)
library(webshot)
library(png)
# conversion to grobs -----------------------------------------------------
htmlwidgets::saveWidget(
  sankey_plot, 
  "temp_sankey.html",
  selfcontained = TRUE
)
webshot(
  "temp_sankey.html", 
  "temp_sankey.png",
  selector = ".sankeyNetwork",
  zoom = 3,
  delay = 0.2,
  vwidth = 500,
  vheight = 600
)
png_data <- readPNG("temp_sankey.png", native = TRUE)
sankey_plot_grob <- rasterGrob(
  png_data,
  interpolate = TRUE,
  width = unit(0.8, "npc"),
  height = unit(1, "npc")
)
file.remove(c("temp_sankey.html", "temp_sankey.png"))
##########################
## polar plot
png(filename = "temp1.png", width = 400, height = 400)
create_polar_plot()
dev.off()
temp_png <- readPNG("temp1.png")
polar_plot_grob <- rasterGrob(temp_png)
file.remove("temp1.png")
class(polar_plot_grob)
# inset -------------------------------------------------------------------
# panel -------------------------------------------------------------------
lay <- rbind(
  c(1, 1, 1, 2, 2),
  c(1, 1, 1, 2, 2),  
  c(5, 5, 5, 2, 2),
  c(3, 3, 3 ,4, 4),
  c(3, 3, 3, 4, 4)
)
plot_list <- list(
  distance_map_plot,
  arranged_plots,
  sankey_plot_grob,
  polar_plot_grob,
  pairwise_dist_plot
)
panel2 <- grid.arrange(
  grobs = plot_list,
  layout_matrix = lay,
  widths =  c(1, 1, 1, 1, 1),
  heights = c(1, 1, 1, 1, 1)
)
# will try saving sepeately for Illustrator -------------------------------
plot_dir <- "plots"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}
ggsave(
  file.path(plot_dir, "distance_map_plot.pdf"),
  plot = distance_map_plot,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(plot_dir, "arranged_plots.pdf"),
  plot = arranged_plots,
  width = 4,
  height = 6,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(plot_dir, "pairwise_dist_plot.pdf"),
  plot = pairwise_dist_plot,
  width = 6,
  height = 4,
  dpi = 300,
  bg = "white"
)
htmlwidgets::saveWidget(
  sankey_plot, 
  "temp_sankey.html",
  selfcontained = TRUE
)
webshot(
  "temp_sankey.html", 
  file.path(plot_dir, "sankey_plot.png"),
  selector = ".sankeyNetwork",
  zoom = 4,
  delay = 0.2,
  vwidth = 1000,
  vheight = 1000
)
png(
  filename = file.path(plot_dir, "polar_plot.png"),
  width = 1000,
  height = 1000,
  res = 300,
  pointsize = 12
)
create_polar_plot()
dev.off()
file.remove("temp_sankey.html")

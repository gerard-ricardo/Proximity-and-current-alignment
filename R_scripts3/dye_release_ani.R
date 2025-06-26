library(magick)
## Add to animation
fps = 2
img_dir <- "E:/1_UQ/Images/2023_palau/2023_04_05_dye_release/drone/centroid"
(image_files <- list.files(path = img_dir, pattern = "*.png", full.names = TRUE))
images <- image_read(image_files)
gif1 <- image_animate(images, fps = fps)
image_write(gif, path = "./plots/compiled_dye1.gif")
 
img_dir <- "E:/1_UQ/Images/2023_palau/2023_04_05_dye_release/drone/subset"
(image_files <- list.files(path = img_dir, pattern = "*.JPG", full.names = TRUE))
images <- image_read(image_files)
gif2 <- image_animate(images, fps = fps)
image_write(gif, path = "./plots/compiled_dye_raw.gif")
# combine panel -----------------------------------------------------------
n_frames <- min(length(gif1), length(gif2))
gif1 <- gif1[1:n_frames]
gif2 <- gif2[1:n_frames]
combined_frames <- image_blank(width = 1, height = 1)
for (i in seq_len(n_frames)) {
  combined_frame <- image_append(c(gif2[i], gif1[i]), stack = FALSE)
  combined_frames <- c(combined_frames, combined_frame)
}
panel_gif <- image_animate(image_join(combined_frames), fps = fps)
image_write(panel_gif, path = "./plots/combined_panel.gif")

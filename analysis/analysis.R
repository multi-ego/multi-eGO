# Title     : TODO
# Objective : TODO
# Created by: Emanuele
# Created on: 26/03/2020

library(tidyverse)
library(reshape2)
#library(hexbin)
#library(gridExtra)
library(grid)
#library(ggforce)

nclust <- read.table(file = "analysis/histo-time.dat", header = FALSE, dec = ".")

  # V1 = timestep
  # V2 = number of monomers
  # V3 = number of dimers
  # V4 = ...
#time <- nclust[,1]
nclust[,1] <- NULL
nclust[,1] <- NULL

s <- 'size'
clust_size <- seq(from = 1, to = (ncol(nclust)))
clust_size <- clust_size + 1
nclust <- t(t(nclust)*clust_size)
clust_size <- paste0(s, "_", clust_size)
colnames(nclust) <- clust_size

# Subset to create an histogram and compare with wet lab results
histogram <- subset(nclust, select = size_2:ncol(nclust))
#is_a_fibril <- subset(nclust, select = size_25:ncol(nclust))
#histogram <- cbind(histogram, total = rowSums(is_a_fibril))
total <- cbind(total = rowSums(histogram))
total_to_plot <- melt(total, varnames = c("time", "fibril"))

# Matrix of the entire clustsize to make the heatmap
matrix_nclust <- melt(nclust, varnames = c("time", "cluster_size"))
colnames(matrix_nclust)[3] <- "cluster_amount_MW"
matrix_nclust[matrix_nclust == 0] <- NA

# Fibril histogram
histo_plot <- ggplot(data = total_to_plot, aes(x = time, y = value)) + geom_col(aes(fill = time)) +
  scale_fill_gradient(low = "#AC80A0", high = "#0471A6", guide = "none", limits = c(275, nrow(total_to_plot))) +
  scale_y_continuous(name = "Fibril MW", breaks = c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500,
                                                    550, 600, 650, 700, 750, 800, 850, 900, 1000)) +
  theme(panel.background = element_rect(fill = "transparent", colour = "black"),
        panel.grid.major = element_line(colour = "grey90"))

# HEAT MAP
heat_clust <- ggplot(data = matrix_nclust, aes(x = time, y = cluster_size)) +
  geom_raster(aes(fill = cluster_amount_MW)) +
  scale_fill_gradient(low = "#AC80A0", high = "#0471A6", na.value = "transparent", guide = "none") +
  scale_y_discrete(breaks = c("size_100", "size_200", "size_300", "size_400", "size_500", "size_600"),
                   name = "Molecules in a Cluster", labels = c("size_100" = "100", "size_200" = "200",
                                                               "size_300" = "300", "size_400" = "400",
                                                               "size_500" = "500", "size_600" = "600")) +
  labs(title = "Heat map of fibrils elongation", fill = "Amount") +
  scale_x_continuous(name = "Time (ns)") +
  theme(panel.background = element_rect(fill = "transparent", colour = "black"),
        panel.grid.major = element_line(colour = "grey90"))

# Merged ones
histo_merge <- ggplot(data = total_to_plot, aes(x = time, y = value)) + geom_col(aes(fill = time), width = 1) +
  scale_fill_gradient(low = "#AC80A0", high = "#0471A6", guide = "none", limits = c(275, nrow(total_to_plot))) +
  scale_y_continuous(name = "Fibril MW") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = "black"),
        panel.grid.major = element_line(colour = "grey90"))

heat_merge <- ggplot(data = matrix_nclust, aes(x = time, y = cluster_size)) +
  geom_raster(aes(fill = cluster_amount_MW)) +
  #geom_raster(aes(fill = cluster_amount_MW)) +
  scale_fill_gradient(low = "#AC80A0", high = "#0471A6", na.value = "transparent", guide = "none") +
  scale_y_discrete(breaks = c("size_100", "size_200", "size_300", "size_400", "size_500", "size_600"),
                   name = "Molecules nmr x MW", labels = c("size_100" = "100", "size_200" = "200",
                                                           "size_300" = "300", "size_400" = "400",
                                                           "size_500" = "500", "size_600" = "600")) +
  scale_x_continuous(name = "Time (ns)") +
  theme(panel.background = element_rect(fill = "transparent", colour = "black"),
        panel.grid.major = element_line(colour = "grey90"))


#heat_clust
png("pd.png")
merge_plot <- grid.newpage()
merge_plot <- grid.draw(rbind(ggplotGrob(histo_merge), ggplotGrob(heat_merge), size = "last"))
dev.off()
merge_plot <- grid.newpage()
merge_plot <- grid.draw(rbind(ggplotGrob(histo_merge), ggplotGrob(heat_merge), size = "last"))
merge_plot
ggsave(plot = merge_plot,"mymerge.png")
heat_zoom <- heat_clust + facet_grid() + coord_cartesian(xlim = c(200,500), ylim = c(0,50)) +
  scale_y_discrete(breaks = c("size_5", "size_10", "size_15", "size_20", "size_25", "size_30", "size_35", "size_40",
                              "size_45", "size_50", name = "none")) +
  scale_x_continuous(breaks = c(200, 250, 275, 300, 350, 400, 450, 500))
merge_zoom <- grid.newpage()
merge_zoom <- grid.draw(rbind(ggplotGrob(heat_merge), ggplotGrob(heat_zoom), size = "last"))
#histo_plot

geom_point(aes(colour = "time"), alpha = 1/50)












# CIMITERO DEGLI ELEFANTI
# QUESTA è UN'IDEA MOLTO CARINA MA STO AVENDO UN ATTIMO DI PROBLEMI
#point_map <- ggplot(data = matrix_nclust, aes(x = time, y = cluster_size)) + geom_point(aes(size = cluster_amount_MW,
#                                                                               colour = cluster_amount_MW),
#                                                                           alpha = 1/50) +
#  scale_y_discrete(breaks = c("size_100", "size_200", "size_300", "size_400", "size_500", "size_600", "size_700"),
#                   name = "Molecules in a Cluster") +
#  scale_radius(guide = "none") +
#  scale_colour_gradient(low = "#AC80A0", high = "#0471A6", na.value = "white", name = "Amount", guide = "none") +
#  theme_bw() +
#  expand_limits(y = 700)




# MERGE
#merge <- ggplot(data = matrix_nclust, aes(x = time, y = cluster_size)) + geom_raster(aes(fill = cluster_amount_MW)) +
#  scale_fill_gradient(low = "#AC80A0", high = "#0471A6", na.value = "white", guide = "none") + # Questo va bene
#  labs(title = "Heat map of fibrils elongation") +
#  geom_point(aes(size = cluster_amount_MW, colour = cluster_amount_MW), alpha = 1/40) +
#  scale_radius(guide = "none") +
#  scale_y_discrete(breaks = c("size_100", "size_200", "size_300", "size_400", "size_500", "size_600"),
#                   name = "Molecules in a Cluster") +
#  scale_colour_gradient(low = "#AC80A0", high = "#0471A6", na.value = "white", guide = "none") +
#  expand_limits(y = 700) +
#  theme_bw()
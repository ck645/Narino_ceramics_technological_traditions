# Disentangling technological traditions: comparative chaines opératoires of 
# painted pre-Hispanic ceramics from Narino, Colombia
# C. Klesner
# 2024



## Broken into relevant sections 
# 1. GMM
# 2. Compositional data - pXRF
# 3. Compositional data - SEM
# 4. Radiocarbon



## Ensure you have installed necessary R packages

# Ensure necessary R packages are installed
packages <- c("dplyr", "rio", "tidyverse", "plyr",
              "Momocs", "stats", "corrr", "rcarbon",
              "ggplot2", "ggbeeswarm", "ggdist", "gghalves", 
              "ggbiplot", "dendextend", "ggtern", "RColorBrewer")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}


## Color conventions for the script - based on brewer palette, Dark2 
# "Capulí" = "#1B9E77", "Piartal" = "#7570B3", "Tuza" = "#D95F02", 
# "Tuza - Red Slip" = "#E7298A", "Piartal - 1" = "#E6AB02", 
# "Piartal - 2" = "#A6761D", #  "Undecorated" = "#666666"


#--------------------------------**GMM**------------------------------------

# The GMM Workflow adapted from structure provided in the workshop: "Geometric 
# Morphometrics and Archaeological Science" by "Dr. Christian Steven Hoggard 
# (University of Southampton, United Kingdom)"


## Load relevant packages for all GMM analysis

library(Momocs)
library(rio)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggbeeswarm)
library(ggdist)
library(gghalves)
library(vegan)
library(dplyr)

## **Creating and Importing GMM Data**

# Note - All data has been digitised using tspUtil and tpsDig2 

tpsdata <- Momocs::import_tps("intervesselcomp.TPS")
GMM_database <- rio::import("GMM_database.csv")


# assigning factors

# Note - V1 = Sample, V3 = Ware, V5 = Site


GMM_database$V1 <- as.factor(GMM_database$V1)
is.factor(GMM_database$V1)

GMM_database$V3 <- as.factor(GMM_database$V3)
is.factor(GMM_database$V3)
summary(GMM_database$V3)

GMM_database$V5 <- as.factor(GMM_database$V5)
is.factor(GMM_database$V5)


## **Creation of the "Out" object**  

# WARNING - you must first ensure that the *tpsdata$coo* names 
# and the database IDs are in the same order

table(names(tpsdata$coo)==GMM_database$V1)
shape <- Out(tpsdata$coo, fac = GMM_database)


# check outlines

Momocs::mosaic(shape, f=GMM_database$V3, pal = pal_qual_Dark2)


## **Outline normalisation**  

# scale, align, rotate, and centre specimens

shapenorm <- coo_center(shape)
shapenorm <- coo_scale(shapenorm)
shapenorm <- coo_close(shapenorm)
shapenorm2 <- shapenorm %>% coo_slidedirection("right") %>% coo_untiltx()


# check outlines:

stack(shapenorm, title = "Stack: Normalised Outlines")
stack(shapenorm2, title = "Stack: Normalised Outlines with coo_slidedirection(right)")


## **Elliptic Fourier Analysis**

calibrate_harmonicpower_efourier(shapenorm2, nb.h = 30)
calibrate_reconstructions_efourier(shapenorm2, range = 1:30)
calibrate_deviations_efourier(shapenorm2)

# generate EFA outlines with harmonics = 99.9% harmonic power

efashape <- efourier(shapenorm2, nb.h = 30, smooth.it = 0, norm = TRUE)


## **Principal Component Analysis**

# use efashape for pca

pcashape <- PCA(efashape)


# determine PC contributions

Scree <- scree(pcashape)


# plot Scree - pick number of PCs that amount to >90%

PC_Scree <- scree_plot(pcashape, nax = 1:5)
PC_Scree + ggtitle("PC Contribution")


# plot PC contribution - pick number of PCs that amount to >90%

PC_Contribution <- PCcontrib(pcashape, nax = 1:3, 
                             sd.r = c(-1, -0.5, 0, 0.5, 1))


# Save PNG file
png(filename = "Figure4.png", width = 2400, height = 1000, res=300)
PCcontrib(pcashape, nax = 1:3, 
          sd.r = c(-1, -0.5, 0, 0.5, 1))
dev.off()


# plot PCA by Ware type - no ellipses

plot_PCA(pcashape, axes = c(1,2), morphospace_position = "range", 
         center_origin = FALSE, zoom = 0.75, chull = FALSE, palette = pal_qual_Paired, 
         title = "PC 1 vs PC 2")

png(filename = "Figure5a.png", width = 2400, height = 1600, res=300)
plot_PCA(pcashape, axes = c(1,2), morphospace_position = "xy", 
         center_origin = FALSE, zoom = 0.85, chull = FALSE, palette = pal_qual_Paired, 
         title = "PC 1 vs PC 2")
dev.off()


#plot PCA by Ware with confidence ellipses = 90%

plot_PCA(pcashape, axes = c(1,2), GMM_database$V3, morphospace_position = "range_axes", 
         center_origin = FALSE, zoom = 0.75, chull = FALSE, palette = pal_qual_Paired, 
         title = "PC 1 vs PC 2") %>% layer_ellipses(conf = 0.9)

plot_PCA(pcashape, axes = c(1,3), GMM_database$V3, morphospace_position = "range_axes", 
         center_origin = FALSE, zoom = 0.9, chull = FALSE, palette = pal_qual_Paired, 
         title = "PC 1 vs PC 3") %>% layer_ellipses(conf = 0.9)

png(filename = "Figure5b.png", width = 2400, height = 1600, res=300)
plot_PCA(pcashape, axes = c(1,2), GMM_database$V3, morphospace_position = "range_axes", 
         center_origin = FALSE, zoom = 0.75, chull = FALSE, 
         pal_manual(c("#1B9E77", "#666666", "#7570B3", "#D95F02", "#E7298A")),  
         title = "PC 1 vs PC 2") %>% layer_ellipses(conf = 0.9)
dev.off()

png(filename = "PC1vs3supplemental.png", width = 2400, height = 1600, res=300)
plot_PCA(pcashape, axes = c(1,3), GMM_database$V3, morphospace_position = "range_axes", 
         center_origin = FALSE, zoom = 0.9, chull = FALSE, 
         pal_manual(c("#1B9E77", "#666666", "#7570B3", "#D95F02", "#E7298A")),  
         title = "PC 1 vs PC 3") %>% layer_ellipses(conf = 0.9)
dev.off()


#boxplot of PC 

boxplot <- boxplot(pcashape, GMM_database$V3, nax = 1:3)
boxplot + scale_fill_brewer(palette = "Paired") + ggtitle("PC contribution sections")


#export PCA scores

file.create('all_scores.txt')
Momocs::export(pcashape, 'all_scores.txt')


## **Discriminant Analysis (LDA/DA/CVA)**  

# First combine all categories into the three wares

GMM_database_filtered <-  GMM_database %>%
  mutate(V3 = case_when(
    V3 == "Capulí - no paint" ~ "Capulí",
    V3 == "Tuza - Red slip" ~ "Tuza",
    TRUE ~ as.character(V3)  # Keep all other values as they are
  ))
GMM_database_filtered$V3 <- as.factor(GMM_database_filtered$V3)

summary(GMM_database_filtered$V3)

dashape99 <- LDA(pcashape, GMM_database_filtered$V3, retain = 0.99)

plot_LDA(dashape99, axes = c(1,2), zoom = 1, chull = FALSE, 
         morphospace_position = "range_axes", 
         pal_manual(c("#1B9E77", "#7570B3", "#D95F02")))

png(filename = "supplemental_LDA.png", width = 2400, height = 1600, res=300)
plot_LDA(dashape99, axes = c(1,2), zoom = 1, chull = FALSE, 
         morphospace_position = "range_axes", 
         pal_manual(c("#1B9E77", "#7570B3", "#D95F02")))
dev.off()


## **Multivariate Analysis of Variance (MANOVA)**

# Interested to see if there is any significant shape difference by wares? 
# Looking at the PC retaining 90% of total variance and 99% of total variance

efashape %>% MANOVA(GMM_database_filtered$V3)
pcashape %>% MANOVA(GMM_database_filtered$V3, test = c("Pillai"), retain = 0.95)
pcashape %>% MANOVA(GMM_database_filtered$V3, test = c("Pillai"), retain = 0.99)


# Which wares differ?

MANOVA_PW(pcashape, GMM_database$V3, retain = 0.99)
pcashape %>% MANOVA_PW(GMM_database_filtered$V3, retain = 0.99)


## **Hierarchical and K-Means Cluster Analysis** 

CLUST(pcashape, retain = 0.90, GMM_database$V3, dist_method = "euclidean", 
      hclust_method = "average", 
      palette = pal_manual(c("#1B9E77", "#666666", "#7570B3", "#D95F02", "#E7298A")))

png(filename = "supplemental_cluster.png", width = 2400, height = 2400, res=300)
CLUST(pcashape, retain = 0.90, GMM_database$V3, dist_method = "euclidean", 
      hclust_method = "average", 
      palette = pal_manual(c("#1B9E77", "#666666", "#7570B3", "#D95F02", "#E7298A")))
dev.off()


## **Mean shapes based on wares**

meanshape <- Out(tpsdata$coo, fac = GMM_database_filtered$V3)

meanshapenorm <- coo_center(meanshape)
meanshapenorm <- coo_scale(meanshapenorm)
meanshapenorm <- coo_close(meanshapenorm)
meanshapenorm2 <- meanshapenorm %>% coo_slidedirection("right") %>% coo_untiltx()

efameanshape <- efourier(meanshapenorm2, nb.h = 30, smooth.it = 0, norm = TRUE)

Meanshapes <- MSHAPES(efameanshape, GMM_database_filtered$V3)

plot_MSHAPES(Meanshapes, size = 0.9)

Out(Meanshapes$shp) %>% panel(names=TRUE)

Out(Meanshapes$shp) %>% panel(names=FALSE)

png(filename = "Figure3d.png", width = 2400, height = 2400, res=300)
Out(Meanshapes$shp) %>% panel(names=FALSE)
dev.off()



## Variance of PCA scores**  

#Import your PCA scores from csv. Csv generated from exported txt file (above) a step needed given the structure of momocs pcashape)

PC_85 <- import('PC_scores.csv')

head(PC_85)

PC_85$V3 <- as.factor(PC_85$V3)

is.factor(PC_85$V3)

PC_85 <- as.data.frame(PC_85)

is.data.frame(PC_85)

#Clean data to only include large ware groups

PC_85 <-  PC_85 %>%
  mutate(V3 = case_when(
    V3 == "Capulí - no paint" ~ "Capulí",
    V3 == "Tuza - Red slip" ~ "Tuza",
    TRUE ~ as.character(V3)  # Keep all other values as they are
  ))
PC_85$V3 <- as.factor(PC_85$V3)

# box and whisker

manual_colors <- c("#1B9E77", "#7570B3", "#D95F02")

PC1 <- ggplot(PC_85, aes(x = V3, y = PC1, fill = V3, color = V3)) +
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .4, 
    .width = 0, 
    justification = -.2, 
    alpha = 0.5,
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA,
    alpha = 0.2,
    color = "black"
  ) +
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .4, 
    alpha = .5,
  ) +
  scale_fill_manual(values = manual_colors) + 
  scale_color_manual(values = manual_colors) +  
  coord_cartesian(xlim = c(1.2, 2.9), clip = "off") 
plot(PC1)

png(filename = "Figure6a.png", width = 2400, height = 1000, res=300)
plot(PC1)
dev.off()


PC2 <- ggplot(PC_85, aes(x = V3, y = PC2, fill = V3, color = V3)) +
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .4, 
    .width = 0, 
    justification = -.2, 
    alpha = 0.5,
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA,
    alpha = 0.2,
    color = "black"
  ) +
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .4, 
    alpha = .5,
  ) +
  scale_fill_manual(values = manual_colors) + 
  scale_color_manual(values = manual_colors) +  
  coord_cartesian(xlim = c(1.2, 2.9), clip = "off") 
plot(PC2)

png(filename = "Figure6b.png", width = 2400, height = 1000, res=300)
plot(PC2)
dev.off()


PC3 <- ggplot(PC_85, aes(x = V3, y = PC3, fill = V3, color = V3)) +
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .4, 
    .width = 0, 
    justification = -.2, 
    alpha = 0.5,
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA,
    alpha = 0.2,
    color = "black"
  ) +
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .4, 
    alpha = .5,
  ) +
  scale_fill_manual(values = manual_colors) + 
  scale_color_manual(values = manual_colors) +  
  coord_cartesian(xlim = c(1.2, 2.9), clip = "off") 
plot(PC3)

png(filename = "Figure6c.png", width = 2400, height = 1000, res=300)
plot(PC3)
dev.off()



## **Identification of Two Piartal groups**

#follows the protocal detailed above, but Piartal vessels are coded as Piartal - 1 or Piartal - 2 based on their PC1 score 

GMM_database_p <- rio::import("GMM_database_piartal.csv")

GMM_database_p$V1 <- as.factor(GMM_database_p$V1)
is.factor(GMM_database_p$V1)

GMM_database_p$V3 <- as.factor(GMM_database_p$V3)
is.factor(GMM_database_p$V3)

GMM_database_p$V5 <- as.factor(GMM_database_p$V5)
is.factor(GMM_database_p$V5)

table(names(tpsdata$coo)==GMM_database_p$V1)

GMM_database_p <-  GMM_database_p %>%
  mutate(V3 = case_when(
    V3 == "Capulí - no paint" ~ "Capulí",
    V3 == "Tuza - Red slip" ~ "Tuza",
    TRUE ~ as.character(V3)  # Keep all other values as they are
  ))
GMM_database_p$V3 <- as.factor(GMM_database_p$V3)

summary(GMM_database_p$V3)

meanshape_p <- Out(tpsdata$coo, fac = GMM_database_p$V3)

meanshapenorm_p <- coo_center(meanshape_p)
meanshapenorm_p <- coo_scale(meanshapenorm_p)
meanshapenorm_p <- coo_close(meanshapenorm_p)
meanshapenorm2_p <- meanshapenorm_p %>% coo_slidedirection("right") %>% coo_untiltx()

efameanshape_p <- efourier(meanshapenorm2_p, nb.h = 30, smooth.it = 0, norm = TRUE)

Meanshapes_p <- MSHAPES(efameanshape_p, GMM_database_p$V3)

plot_MSHAPES(Meanshapes_p, size = 0.9, palette= pal_manual(c("blue", "red")))

png(filename = "Figure7.png", width = 3600, height = 3600, res=300)
plot_MSHAPES(Meanshapes_p, size = 0.9, palette= pal_manual(c("blue", "red")))
dev.off()

Out(Meanshapes_p$shp) %>% panel(names=TRUE)

plot_PCA(pcashape, axes = c(1,2), GMM_database_p$V3, morphospace_position = "range_axes", 
         center_origin = FALSE, zoom = 0.9, chull = FALSE, palette = pal_qual_Paired, 
         title = "PC 1 vs PC 2") %>% layer_ellipses(conf = 0.9) 

plot_PCA(pcashape, axes = c(1,3), GMM_database_p$V3, morphospace_position = "range_axes", 
         center_origin = FALSE, zoom = 0.9, chull = FALSE, palette = pal_qual_Paired, 
         title = "PC 1 vs PC 3") %>% layer_ellipses(conf = 0.9)

MANOVA_PW(pcashape, GMM_database_p$V3, retain = 0.90)

pcashape %>% MANOVA_PW(GMM_database_p$V3, retain = 0.99)



# boxplots of metric measurements

metric_measurements <- import('metric_vessel_measurements.csv')

head(metric_measurements)
metric_measurements$Ware <- as.factor(metric_measurements$Ware)
is.factor(metric_measurements$Ware)

metric_measurements <- as.data.frame(metric_measurements)
is.data.frame(metric_measurements)

metric_measurements <-  metric_measurements %>%
  mutate(Ware = case_when(
    Ware == "Capulí - no paint" ~ "Capulí",
    Ware == "Capuli" ~ "Capulí",
    Ware == "Tuza - Red slip" ~ "Tuza",
    Ware == "Piartal - 1" ~ "Piartal",
    Ware == "Piartal - 2" ~ "Piartal",
    TRUE ~ as.character(Ware)  # Keep all other values as they are
  ))
metric_measurements$Ware <- as.factor(metric_measurements$Ware)

is.numeric(metric_measurements$`Height (cm)`)

Height <- ggplot(metric_measurements, aes(x = Ware, y = metric_measurements$`Height (cm)`, 
                                          fill = Ware, color = Ware)) +
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .4, 
    .width = 0, 
    justification = -.2, 
    alpha = 0.5,
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA,
    alpha = 0.2,
    color = "black"
  ) +
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .4, 
    alpha = .5,
  ) +
  scale_fill_manual(values = manual_colors) + 
  scale_color_manual(values = manual_colors) +  
  coord_cartesian(xlim = c(1.2, 2.9), clip = "off") 
plot(Height)

png(filename = "Figure3a.png", width = 2400, height = 1000, res=300)
plot(Height)
dev.off()


Base <- ggplot(metric_measurements, aes(x = Ware, 
                                        y = metric_measurements$`Base Radius (cm)`, 
                                        fill = Ware, color = Ware)) +
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .4, 
    .width = 0, 
    justification = -.2, 
    alpha = 0.5,
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA,
    alpha = 0.2,
    color = "black"
  ) +
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .4, 
    alpha = .5,
  ) +
  scale_fill_manual(values = manual_colors) + 
  scale_color_manual(values = manual_colors) +  
  coord_cartesian(xlim = c(1.2, 2.9), clip = "off") 
plot(Base)

png(filename = "Figure3b.png", width = 2400, height = 1000, res=300)
plot(Base)
dev.off()


Rim <- ggplot(metric_measurements, aes(x = Ware, y = `Rim Radius (cm)`, 
                                       fill = Ware, color = Ware)) +
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .4, 
    .width = 0, 
    justification = -.2, 
    alpha = 0.5,
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA,
    alpha = 0.2,
    color = "black"
  ) +
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .4, 
    alpha = .5,
  ) +
  scale_fill_manual(values = manual_colors) + 
  scale_color_manual(values = manual_colors) +  
  coord_cartesian(xlim = c(1.2, 2.9), clip = "off") 
plot(Rim)

png(filename = "Figure3c.png", width = 2400, height = 1000, res=300)
plot(Rim)
dev.off()




#----------------------**Compositional data - pXRF**---------------------------

library(rio)
library(dplyr)
library(ggplot2)
library(ggbiplot)
library(RColorBrewer)
library(stats)
library(dendextend)

# Importing the dataset
database <- rio::import("compdata_piartal.csv")


## **Cleaning the dataset**

database$Sample <- as.factor(database$Sample)
database$Collection <- as.factor(database$Collection)
database$Ware <- as.factor(database$Ware)
database$tomb <- as.factor(database$tomb)
database$component <- as.factor(database$component)
database$Mg <- as.numeric(database$Mg)
database$Ni <- as.numeric(database$Ni)
database$Ti <- as.numeric(database$Ti)
database$Nb <- as.numeric(database$Nb)
database$Pb <- as.numeric(database$Pb)


# remove outliers/non-authentic pieces
database_filtered <- database %>%
  filter(Sample != "CA230346", Sample != "CA230351", Sample != "CA230298_1", Sample != "CA230298_2")

database_filtered <-  database_filtered %>%
  mutate(Collection = case_when(
    Collection == "Bomboná" ~ "Various - Bomboná",
    Collection == "MO" ~ "Various - Museo del Oro",
    Collection == "Union del Sur - Site 8" ~ "Catambuco - Site 8",
    Collection  == "San Damián" ~ "Catambuco - San Damián",
    Collection == "Union del Sur - Site 137" ~ "El Porvenir - Site 137",
    TRUE ~ as.character(Collection)  # Keep all other values as they are
  ))
database_filtered$Collection <- as.factor(database_filtered$Collection)
summary(database_filtered$Collection)


# Substituting 2/3 value of LLOD for elements not routinely detected (Mg, Ni, Ti, Nb, Pb)

detection_limit_Mg <- 9769
replacement_value_Mg <- (2/3) * detection_limit_Mg

detection_limit_Ni <- 15
replacement_value_Ni <- (2/3) * detection_limit_Ni

detection_limit_Ti <- 3031
replacement_value_Ti <- (2/3) * detection_limit_Ti

detection_limit_Nb <- 8
replacement_value_Nb <- (2/3) * detection_limit_Nb

detection_limit_Pb <- 15
replacement_value_Pb <- (2/3) * detection_limit_Pb

database_filtered <- database_filtered %>%
  mutate(Mg = ifelse(is.na(Mg), replacement_value_Mg, Mg))
database_filtered <- database_filtered %>%
  mutate(Ni = ifelse(is.na(Ni), replacement_value_Ni, Ni))
database_filtered <- database_filtered %>%
  mutate(Ti = ifelse(is.na(Ti), replacement_value_Ti, Ti))
database_filtered <- database_filtered %>%
  mutate(Nb = ifelse(is.na(Nb), replacement_value_Nb, Nb))
database_filtered <- database_filtered %>%
  mutate(Pb = ifelse(is.na(Pb), replacement_value_Pb, Pb))

database_filtered_p <- database_filtered


# cleaning ware groups and collections
database_filtered <-  database_filtered %>%
  mutate(Ware = case_when(
    Ware == "Capuli" ~ "Capulí",
    Ware == "Piartal?" ~ "Piartal",
    Ware == "Piartal - 1" ~ "Piartal",
    Ware == "Piartal - 2" ~ "Piartal",
    Ware == "Tuza?" ~ "Tuza",
    Ware == "undec." ~ "Undecorated",
    Ware == "undecorated" ~ "Undecorated",
    Ware == "undecorated?" ~ "Undecorated",
    Ware == "red slip" ~ "Tuza - Red Slip",
    TRUE ~ as.character(Ware)  # Keep all other values as they are
  ))
database_filtered$Ware <- as.factor(database_filtered$Ware)
summary(database_filtered$Ware)


# Subsetting the dataframe to isolate the fabric component
fabric_dataset <- filter(database_filtered, component %in% c("fabric"))


## **plotting compositional biplots**

# select the wares to plot with 90% confidence ellipses
ellipse_data_lim <- fabric_dataset %>% 
  filter(Ware %in% c("Capulí", "Piartal", "Tuza"))


Fe2O3_vs_CaO_ellipse_lim <- ggplot(fabric_dataset, aes(x = Fe2O3, y = CaO, color = Ware, shape = Collection)) +
  geom_point(size=3) +
  # Ellipse outline, allow color to be determined by 'color = Ware' mapping
  stat_ellipse(data = ellipse_data_lim, aes(group = Ware, color = Ware), level = 0.90, geom = "path") +
  theme_minimal() +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Fe2O3 (%)", y = "CaO (%)") +
  scale_color_manual(name = "Ware", values = c("Capulí" = "#1B9E77", "Piartal" = "#7570B3", 
              "Tuza" = "#D95F02", "Tuza - Red Slip" = "#E7298A", "Undecorated" = "#666666"))+ 
  scale_shape_discrete(name = "Site") +
  theme(legend.position = "right")

Fe2O3_vs_CaO_ellipse_lim + theme(axis.line = element_line(color = "black", linewidth = 1, linetype = 1)) 

# Save PNG file
png(filename = "Figure8.png", width = 2400, height = 1600, res=300)
plot(Fe2O3_vs_CaO_ellipse_lim + theme(axis.line = element_line(color = "black", linewidth = 1, linetype = 1)))
dev.off()


## ** Compositional PCA**

# Select relevant elements
pc <- prcomp(fabric_dataset[,c("Al", "K", "Ca", "Ti", "Mn", "Fe", "Zn", "Rb", "Y", "Zr", "Nb", "Pb")],
             center = TRUE,
             scale. = TRUE) 

attributes(pc) 
pc.dataframe <- as.data.frame(pc$x)
pc.dataframe <- data.frame(pc.dataframe, Ware = fabric_dataset$Ware, Site = fabric_dataset$Collection, Sample = fabric_dataset$Sample)
pc.dataframe <-  pc.dataframe %>%
  mutate(Ware = case_when(
    Ware == "Tuza - Red Slip" ~ "Tuza",
    TRUE ~ as.character(Ware)  # Keep all other values as they are
  ))
pc.dataframe$Ware <- as.factor(pc.dataframe$Ware)

# compute the variance explained
explained_variance <- pc$sdev^2
proportion_variance <- explained_variance / sum(explained_variance)
list(proportion_variance)


# create the loadings for the element vectors, and scale them for visualization
loadings <- as.data.frame(pc$rotation)
loadings_scaled <- loadings*8
loadings_scaled$Variable <- rownames(loadings)

# select the wares to plot with 90% confidence ellipses
ellipse_pc.dataframe <- pc.dataframe %>% 
  filter(Ware %in% c("Capulí", "Piartal", "Tuza")) 


g2 <- ggplot(pc.dataframe, aes(x = PC1, y = PC2, color = Ware)) +
  geom_point(size=2) +
  theme_minimal() +
  stat_ellipse(data = ellipse_pc.dataframe, aes(group = Ware, color = Ware), level = 0.90, geom = "path") +
  geom_segment(data = loadings_scaled, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text(data = loadings_scaled, aes(x = PC1, y = PC2, label = Variable), hjust = 0.3, vjust = 1, size = 5, color = "black") +
  labs(title = "PCA Biplot for select elements - 90% confidence ellipses", x = "PC1 (23.1%)", y = "PC2 (14.7%)") +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  scale_color_manual(name = "Ware", values = c("Capulí" = "#1B9E77", "Piartal" = "#7570B3", "Tuza" = "#D95F02", "Undecorated" = "#666666"))+ 
  scale_shape_discrete(name = "Site") +
  theme(legend.position = "right")

print(g2) + theme(axis.line = element_line(color = "black", linewidth = 1, linetype = 1))

# Save PNG file
png(filename = "Figure9.png", width = 2400, height = 1600, res=300)
plot((g2) + theme(axis.line = element_line(color = "black", linewidth = 1, linetype = 1)))
dev.off()


## Examine compositional behavior based on the two different Piartal groups

summary(database_filtered_p$Ware)

# Subsetting the dataframe to look at the fabric
fabric_p_dataset <- filter(database_filtered_p, component %in% c("fabric"))

# Subset further to isolate just the decorated samples
fabric_p_dataset <- filter(fabric_p_dataset, Ware %in% c("Capulí", "Piartal", "Piartal - 1", "Piartal - 2", "Tuza", "Tuza - Red Slip"))


#plotting PC by Piartal groups (following same protocol as above)

pc_p <- prcomp(fabric_p_dataset[,c("Al", "K", "Ca", "Ti", "Mn", "Fe", "Zn", "Rb", "Y", "Zr", "Nb", "Pb")],
               center = TRUE,
               scale. = TRUE) 

attributes(pc_p) 

pc_p.dataframe <- as.data.frame(pc_p$x)
pc_p.dataframe <- data.frame(pc_p.dataframe, Ware = fabric_p_dataset$Ware, Site = fabric_p_dataset$Collection, Sample = fabric_p_dataset$Sample)

summary(pc_p.dataframe$Ware)

pc_p.dataframe <-  pc_p.dataframe %>%
  mutate(Ware = case_when(
    Ware == "Tuza - Red Slip" ~ "Tuza",
    TRUE ~ as.character(Ware)  # Keep all other values as they are
  ))
pc_p.dataframe$Ware <- as.factor(pc_p.dataframe$Ware)
summary(pc_p.dataframe$Ware)

summary(pc_p.dataframe$Site)

# Loadings scaled for visualization
loadings_p <- as.data.frame(pc_p$rotation)
loadings_p_scaled <- loadings_p*8
loadings_p_scaled$Variable <- rownames(loadings_p)

ellipse_pc_p.dataframe <- pc_p.dataframe %>% 
  filter(Ware %in% c("Piartal - 1", "Piartal", "Piartal - 2")) 

g2_p <- ggplot(pc_p.dataframe, aes(x = PC1, y = PC2, color = Ware)) +
  geom_point(size=2) +
  theme_minimal() +
  stat_ellipse(data = ellipse_pc_p.dataframe, aes(group = Ware, color = Ware), level = 0.90, geom = "path") +
  geom_segment(data = loadings_scaled, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text(data = loadings_scaled, aes(x = PC1, y = PC2, label = Variable), hjust = 0.3, vjust = 1, size = 5, color = "black") +
  labs(title = "PCA Biplot for select elements - 90% confidence ellipses", x = "PC1", y = "PC2") +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  scale_color_brewer(name = "Ware", palette = "Set1") +
  scale_shape_discrete(name = "Site") +
  theme(legend.position = "right")

print(g2_p) + theme(axis.line = element_line(color = "black", linewidth = 1, linetype = 1))

#Plot PCA - biplot isolating Piartal

g4_e <- ggplot(pc_p.dataframe, aes(x = PC1, y = PC2, color = Ware)) +
  geom_point(size=3) +
  theme_minimal() +
  stat_ellipse(data = ellipse_pc_p.dataframe, aes(group = Ware, color = Ware), level = 0.90, geom = "path") +
  geom_segment(data = loadings_scaled, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text(data = loadings_scaled, aes(x = PC1, y = PC2, label = Variable), hjust = 0.3, vjust = 1, size = 5, color = "black") +
  labs(title = "PCA Biplot for select elements - Piartal - 90% confidence ellipses", x = "PC1", y = "PC2") +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  scale_shape_discrete(name = "Ware") +
  scale_color_manual(
    name = "Ware",
    values = c("Piartal" = "#7570B3", "Piartal - 1" = "#E6AB02", "Piartal - 2" = "#A6761D", "Tuza" = "darkgrey", "Capulí" = "lightgrey"),
    breaks = c("Piartal", "Piartal - 1", "Piartal - 2")
  ) +
  theme(legend.position = "right")

print(g4_e) + theme(axis.line = element_line(color = "black", linewidth = 1, linetype = 1))

# Save PNG file
png(filename = "Figure10a.png", width = 2400, height = 2000, res=300)
plot((g4_e) + theme(axis.line = element_line(color = "black", linewidth = 1, linetype = 1)))
dev.off()

# without loadings
shapes3 <-  c(3, 17, 15, 19, 4)

g4 <- ggplot(pc_p.dataframe, aes(x = PC1, y = PC2, color = Ware, shape = Ware)) +
  geom_point(size=2, linewidth=1) +
  theme_minimal() +
  stat_ellipse(data = ellipse_pc_p.dataframe, aes(group = Ware, color = Ware), level = 0.90, geom = "path") +
  labs(title = "PCA Biplot for select elements - Piartal - 90% confidence ellipses", x = "PC1", y = "PC2") +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  scale_shape_discrete(name = "Ware") +
  scale_shape_manual(values = shapes3, name = "Sample") +
  scale_color_manual(
    name = "Ware",
    values = c("Piartal" = "#7570B3", "Piartal - 1" = "#E6AB02", "Piartal - 2" = "#A6761D", "Tuza" = "darkgrey", "Capulí" = "lightgrey"),
    breaks = c("Piartal", "Piartal - 1", "Piartal - 2")
  ) +
  theme(legend.position = "right")

print(g4) + theme(axis.line = element_line(color = "black", linewidth = 1, linetype = 1))

# Save PNG file
png(filename = "Figure10b.png", width = 2400, height = 2000, res=300)
plot((g4) + theme(axis.line = element_line(color = "black", linewidth = 1, linetype = 1)))
dev.off()


## **Hierarchical Cluster Analysis** 

# Select the columns for clustering
fabric_cluster <- filter(fabric_p_dataset, Ware %in% c("Capulí", "Piartal", "Piartal - 1", 
                                      "Piartal - 2", "Tuza", "Tuza - Red Slip"))
fabric_cluster <- fabric_cluster %>%
  select(Sample, Ware, Al, K, Ca, Ti, Fe, Mn, Zn, Rb, Y, Zr, Nb, Pb)

# Standardize the data 
data_for_clustering_scaled <- as.data.frame(scale(select(fabric_cluster, -Sample, -Ware)))

# Perform hierarchical clustering
distance_matrix <- dist(data_for_clustering_scaled)
hc <- hclust(distance_matrix, method = "complete")

# Cluster analysis by Ware
dend <- as.dendrogram(hc)
labels(dend) <- fabric_cluster$Sample[order.dendrogram(dend)]

ware_colors <- c(
  "Capulí" = "#1B9E77",
  "Piartal" = "#7570B3",
  "Piartal - 1" = "#E6AB02",
  "Piartal - 2" = "#A6761D",
  "Tuza" = "#D95F02",
  "Tuza - Red Slip" = "#E7298A"
)
fabric_cluster$Ware <- factor(fabric_cluster$Ware, levels = names(ware_colors))

labels_colors(dend) <- ware_colors[fabric_cluster$Ware[order.dendrogram(dend)]]

dend <- set(dend, "labels_cex", 0.6)  # Adjust this value as needed
plot(dend, main = "Hierarchical Clustering of Samples")
legend("topright", legend = names(ware_colors), fill = ware_colors, title = "Ware Group", cex = 0.5)

# Save PNG file
png(filename = "supplemental_hc.png", width = 2400, height = 1600, res=300)
plot(dend, main = "Hierarchical Clustering of Samples")
legend("topright", legend = names(ware_colors), fill = ware_colors, title = "Ware Group", cex = 0.5)
dev.off()

#----------------------**Compositional data - SEM**-------------------------------

library(rio)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggtern)
library(plyr)
library(devtools) 
library(ggbiplot) 
library("corrr")


# Importing the dataset
SEM_database <- rio::import("inclusions.csv")
feldspars <- rio::import("Fldspr.csv")


## **Cleaning the dataset**
SEM_database$Sample <- as.factor(SEM_database$Sample)
SEM_database$Ware <- as.factor(SEM_database$Ware)
SEM_database$component <- as.factor(SEM_database$component)

feldspars$Sample <- as.factor(feldspars$Sample)
feldspars$Ware <- as.factor(feldspars$Ware)

# checking groups
summary(SEM_database$Ware)
summary(SEM_database$Sample)
summary(SEM_database$component)

SEM_database <-  SEM_database %>%
  mutate(Ware = case_when(
    Ware == "red slip" ~ "Tuza - Red Slip",
    TRUE ~ as.character(Ware)  # Keep all other values as they are
  ))
SEM_database$Ware <- as.factor(SEM_database$Ware)
summary(SEM_database$Ware)
summary(SEM_database$Sample)

SEM_database <-  SEM_database %>%
  mutate(component = case_when(
    component == "Qtz" ~ "quartz",
    component == "spherical inclusion" ~ "quartz",
    component == "Quartz" ~ "quartz",
    component == "porous inclusion" ~ "Other minerals",
    component == "inclusion trapped in glassy inclusion" ~ "Other minerals",
    component == 'inclusion' ~ 'Other minerals',
    component == 'inclusion with Fe-Ti zoning' ~ 'Fe-Ti inclusion',
    TRUE ~ as.character(component)  # Keep all other values as they are
  ))
SEM_database$component <- as.factor(SEM_database$component)
summary(SEM_database$component)

feldspars <-  feldspars %>%
  mutate(Ware = case_when(
    Ware == "red slip" ~ "Tuza - Red Slip",
    Ware == "Capul\xed" ~ "Capulí",
    TRUE ~ as.character(Ware)  # Keep all other values as they are
  ))
feldspars$Ware <- as.factor(feldspars$Ware)
summary(feldspars$Ware)


## **Subsetting the dataframe** 

# Creating groups based on 'component' factors
fabric_SEM_dataset <- filter(SEM_database, component %in% c("fabric"))
glassy_dataset <- filter(SEM_database, component %in% c("glassy inclusion"))
feldspars <- filter(SEM_database, component %in% c("feldspar"))
FeTi_dataset <- filter(SEM_database, component %in% c("Fe-Ti inclusion"))
generalinclusion_dataset <- filter(SEM_database, component %in% c("inclusion", "feldspar"))
mostinclusion_dataset <- filter(SEM_database, component %in% c("inclusion", "feldspar", "glassy inclusion"))


## **Plotting compositional biplots**

shapes <- c(1, 2, 8, 6, 15, 0, 16, 17)

Si_vs_Al <- ggplot(glassy_dataset, aes(x = SiO2, y = Al2O3, color = Ware, shape = Sample)) +
  geom_point(size=2) +
  theme_minimal() +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  labs(title = "Si vs Al", x = "SiO2 (wt%)", y = "Al2O3 (wt%)") +
  scale_color_manual(name = "Ware", 
                     values = c("Capulí" = "#1B9E77", "Piartal" = "#7570B3", 
                                "Tuza" = "#D95F02", "Tuza - Red Slip" = "#E7298A"))+
  scale_shape_manual(values = shapes, name = "Sample") +
  theme(legend.position = "right")

Si_vs_Al + theme(axis.line = element_line(color = "black", linewidth = 1, linetype = 1))

png(filename = "Figure16.png", width = 2400, height = 1600, res=300)
plot(Si_vs_Al + theme(axis.line = element_line(color = "black", linewidth = 1, linetype = 1)))
dev.off()


## **Plotting ternary diagrams** 

feldsparternary <- ggtern(data=feldspars, aes(x=Ab, y=An, z=Or, fill=Ware)) +
  geom_point(shape = 21, color = "black", size = 2.5, stroke = 0.5, alpha = 0.8) +  
  theme_minimal() +
  theme(axis.line = element_line(color = "black", linewidth = 1, linetype = 1)) +
  labs(title = "Feldspar Ternary", x = "Ab", y = "An", z = "Or") +
  scale_fill_manual(name = "Ware", 
                    values = c("Capulí" = "#1B9E77", "Piartal" = "#7570B3", 
                               "Tuza" = "#D95F02", "Tuza - Red Slip" = "#E7298A"))+
  theme(legend.position = "right") +
  theme(
    tern.axis.line.T = element_line(color = "black", size = 1.2),
    tern.axis.line.L = element_line(color = "black", size = 1.2),
    tern.axis.line.R = element_line(color = "black", size = 1.2),
    tern.panel.grid.major.T = element_blank(),
    tern.panel.grid.major.L = element_blank(),
    tern.panel.grid.major.R = element_blank(),
    tern.panel.grid.minor.T = element_blank(),
    tern.panel.grid.minor.L = element_blank(),
    tern.panel.grid.minor.R = element_blank(), 
    panel.background = element_blank(),  
    plot.background = element_blank()  
  )

plot(feldsparternary)

png(filename = "Figure13.png", width = 2400, height = 2400, res=300)
plot(feldsparternary)
dev.off()


Capulifeldsparternary <- ggtern(data=feldspars, aes(x=Ab, y=An, z=Or, fill=Ware)) +
  geom_point(shape = 21, color = "black", size = 3, stroke = 0.5) +  
  scale_fill_manual(values = c("Capulí" = "#1B9E77", "Piartal" = "transparent", 
                               "Tuza" = "transparent","Tuza - Red Slip" = "transparent")) +
  theme_minimal() +
  theme(axis.line = element_line(color = "black", linewidth = 1, linetype = 1)) +
  labs(title = "Feldspar Ternary", x = "Ab", y = "An", z = "Or") +
  scale_color_brewer(name = "Ware", palette = "Set1") +
  theme(legend.position = "right") +
  theme(
    tern.axis.line.T = element_line(color = "black", size = 1.2),
    tern.axis.line.L = element_line(color = "black", size = 1.2),
    tern.axis.line.R = element_line(color = "black", size = 1.2),
    tern.panel.grid.major.T = element_blank(),
    tern.panel.grid.major.L = element_blank(),
    tern.panel.grid.major.R = element_blank(),
    tern.panel.grid.minor.T = element_blank(),
    tern.panel.grid.minor.L = element_blank(),
    tern.panel.grid.minor.R = element_blank(), 
    panel.background = element_blank(),  
    plot.background = element_blank()  
  )

png(filename = "Figure14a.png", width = 2400, height = 2400, res=300)
plot(Capulifeldsparternary)
dev.off()


Piartalfeldsparternary <- ggtern(data=feldspars, aes(x=Ab, y=An, z=Or, fill=Ware)) +
  geom_point(shape = 21, color = "black", size = 3, stroke = 0.5) +  
  scale_fill_manual(values = c("Capulí" = "transparent", "Piartal" = "#7570B3", 
                               "Tuza" = "transparent","Tuza - Red Slip" = "transparent")) +
  theme_minimal() +
  theme(axis.line = element_line(color = "black", linewidth = 1, linetype = 1)) +
  labs(title = "Feldspar Ternary", x = "Ab", y = "An", z = "Or") +
  scale_color_brewer(name = "Ware", palette = "Set1") +
  theme(legend.position = "right") +
  theme(
    tern.axis.line.T = element_line(color = "black", size = 1.2),
    tern.axis.line.L = element_line(color = "black", size = 1.2),
    tern.axis.line.R = element_line(color = "black", size = 1.2),
    tern.panel.grid.major.T = element_blank(),
    tern.panel.grid.major.L = element_blank(),
    tern.panel.grid.major.R = element_blank(),
    tern.panel.grid.minor.T = element_blank(),
    tern.panel.grid.minor.L = element_blank(),
    tern.panel.grid.minor.R = element_blank(), 
    panel.background = element_blank(),  
    plot.background = element_blank()  
  )

png(filename = "Figure14b.png", width = 2400, height = 2400, res=300)
plot(Piartalfeldsparternary)
dev.off()


Tuzafeldsparternary <- ggtern(data=feldspars, aes(x=Ab, y=An, z=Or, fill=Ware)) +
  geom_point(shape = 21, color = "black", size = 3, stroke = 0.5) +  
  scale_fill_manual(values = c("Capulí" = "transparent", "Piartal" = "transparent", 
                                "Tuza" = "#D95F02","Tuza - Red Slip" = "#E7298A")) +
  theme_minimal() +
  theme(axis.line = element_line(color = "black", linewidth = 1, linetype = 1)) +
  labs(title = "Feldspar Ternary", x = "Ab", y = "An", z = "Or") +
  theme(legend.position = "right") +
  theme(
    tern.axis.line.T = element_line(color = "black", size = 1.2),
    tern.axis.line.L = element_line(color = "black", size = 1.2),
    tern.axis.line.R = element_line(color = "black", size = 1.2),
    tern.panel.grid.major.T = element_blank(),
    tern.panel.grid.major.L = element_blank(),
    tern.panel.grid.major.R = element_blank(),
    tern.panel.grid.minor.T = element_blank(),
    tern.panel.grid.minor.L = element_blank(),
    tern.panel.grid.minor.R = element_blank(), 
    panel.background = element_blank(),  
    plot.background = element_blank()  
  )

plot(Tuzafeldsparternary)

png(filename = "Figure14c.png", width = 2400, height = 2400, res=300)
plot(Tuzafeldsparternary)
dev.off()


## **Create histograms of feldspars anorthite and albite contributions**

#compute mean anorthite (An) contibutions
mean_An <- ddply(feldspars, "Ware", summarise, mean = mean(An, na.rm = TRUE))
head(mean_An)

anorthite <- ggplot(data=feldspars, aes(x=An, fill=Ware)) + 
  geom_histogram(color="black", alpha=0.8,) +
  facet_grid(Ware ~ .) +
  geom_vline(data=mean_An, aes(xintercept=mean),
             linetype="dashed") +
  scale_fill_manual(name = "Ware", 
                    values = c("Capulí" = "#1B9E77", "Piartal" = "#7570B3", 
                               "Tuza" = "#D95F02", "Tuza - Red Slip" = "#E7298A")) +
  labs(title = "Anorthite", x = "Anorthite") +
  theme(panel.background = element_blank(),  
  plot.background = element_blank())

plot(anorthite)

png(filename = "anorthite_supplemental.png", width = 2400, height = 2400, res=300)
plot(anorthite)
dev.off()

#compute mean albite (Ab) contibutions
mean_Ab <- ddply(feldspars, "Ware", summarise, mean = mean(Ab, na.rm = TRUE))
head(mean_Ab)

albite <- ggplot(data=feldspars, aes(x=Ab, fill=Ware)) + 
  geom_histogram(color="black", alpha=0.8,) +
  facet_grid(Ware ~ .) +
  geom_vline(data=mean_Ab, aes(xintercept=mean),
             linetype="dashed") +
  scale_fill_manual(name = "Ware", 
                    values = c("Capulí" = "#1B9E77", "Piartal" = "#7570B3", 
                               "Tuza" = "#D95F02", "Tuza - Red Slip" = "#E7298A")) +
  labs(title = "Albite", x = "Albite") +
  theme(panel.background = element_blank(),  
        plot.background = element_blank())

plot(albite)

png(filename = "albite_supplemental.png", width = 2400, height = 2400, res=300)
plot(albite)
dev.off()


# isolate the feldspars just in the Tuza ware ceramics
feldspars_Tuza <- filter(feldspars, Ware %in% c("Tuza"))

mean_T_Ab <- ddply(feldspars_Tuza, "Sample", summarise, mean = mean(Ab, na.rm = TRUE))
head(mean_T_Ab)

albite_Tuza <- ggplot(data=feldspars_Tuza, aes(x=Ab)) + 
  geom_histogram(color="black", fill="#D95F02", alpha=1,) +
  facet_grid(Sample ~ .) +
  geom_vline(data=mean_T_Ab, aes(xintercept=mean),
             linetype="dashed") +
  labs(title = "Albite - Tuza samples", x = "Albite") +
  theme(panel.background = element_blank(),  
        plot.background = element_blank())

plot(albite_Tuza)

png(filename = "albite_Tuza_supplemental.png", width = 2400, height = 2400, res=300)
plot(albite_Tuza)
dev.off()


## PCA of the inclusions

# select all inclusions but exclude fabric
PCA_dataset <- filter(SEM_database, 
            component %in% c("Other minerals", "feldspar", 
            "glassy inclusion", "Fe-Ti inclusion",  
            "quartz", "rock fragment", "slip inclusion"))

# PCA on all nine major oxides (PCA follows procedures above)
pcall <- prcomp(PCA_dataset[,c("SiO2", "Na2O", "MgO", "Al2O3", "P2O5", "K2O", 
                               "CaO", "TiO2", "FeO")],
                center = TRUE,
                scale. = TRUE) 
attributes(pcall) 

PC_all  <- ggbiplot(pcall,
                    obs.scale = 1,
                    var.scale = 1,
                    groups = PCA_dataset$component,
                    ellipse = TRUE, 
                    ellipse.prob = 0.90) 

PC_all <- PC_all + scale_color_discrete(name = '') 

print(PC_all)

pcall.dataframe <- as.data.frame(pcall$x)

pcall.dataframe <- data.frame(pcall.dataframe, Ware = PCA_dataset$Ware, 
      component = PCA_dataset$component, Sample = PCA_dataset$Sample, 
      Sample_ID = PCA_dataset$`n= (for fabric) or analysis # (for inclusions)`)

summary(pcall.dataframe$component)

file.create('all_scores_SEM.txt')
export(pcall.dataframe, 'all_scores_SEM.txt')

# calculate variance explained
explained_variance_inclusions <- pcall$sdev^2
proportion_variance_inclusions <- explained_variance_inclusions / sum(explained_variance_inclusions)
list(proportion_variance_inclusions)

# Designate loadings for element vectors, scaled for visualization
loadingsall <- as.data.frame(pcall$rotation)
loadingsall_scaled <- loadingsall*7
loadingsall_scaled$Variable <- rownames(loadingsall)

#Plot PCA by ware
PC_all.1 <- ggplot(pcall.dataframe, aes(x = PC1, y = PC2, fill=Ware)) +
  geom_point(shape = 21, color = "black", size = 2.5, stroke = 0.5, alpha = 0.8) +  
  scale_fill_manual(name = "Ware", 
                    values = c("Capulí" = "#1B9E77", "Piartal" = "#7570B3", 
                               "Tuza" = "#D95F02", "Tuza - Red Slip" = "#E7298A"))+
  theme_minimal() +
  labs(title = "PCA Biplot for major elements", x = "PC1 (34.1%)", y = "PC2 (21.3%)") +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  xlim(-6, 4) +
  ylim(-4, 3) +
  theme(legend.position = "right")

print(PC_all.1) + theme(axis.line = element_line(color = "black", linewidth = 1, linetype = 1))

png(filename = "Figure12a.png", width = 3600, height = 2400, res=300)
plot(PC_all.1) + theme(axis.line = element_line(color = "black", linewidth = 1, linetype = 1))
dev.off()


#Plot PCA by inclusion type 
PC_all.2 <- ggplot(pcall.dataframe, aes(x = PC1, y = PC2, color=component)) +
  geom_point(size = 2) +  
  scale_color_brewer(palette = "Set2") +
  theme_minimal() +
  geom_segment(data = loadingsall_scaled, aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text(data = loadingsall_scaled, aes(x = PC1, y = PC2, label = Variable), 
            hjust = 0.3, vjust = 1, size = 5, color = "black") +
  labs(title = "PCA Biplot for major elements", x = "PC1 (34.1%)", y = "PC2 (21.3%)") +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  xlim(-6, 4) +
  ylim(-4, 3) +
  theme(legend.position = "right")

print(PC_all.2) + theme(axis.line = element_line(color = "black", linewidth = 1, linetype = 1))

png(filename = "Figure12b.png", width = 3600, height = 2400, res=300)
plot(PC_all.2) + theme(axis.line = element_line(color = "black", linewidth = 1, linetype = 1))
dev.off()


#----------------------**Radiocarbon**-----------------------------------------

#Load the relevant packages

library(rcarbon)
library(rio)
library(dplyr)

# import the data
alldates <- rio::import("radiocarbon.csv")

## **Clean the dataframe**
head(alldates)

alldates$Municipality <- as.factor(alldates$Municipality)
alldates$`Ceramic Ware identification` <- as.factor(alldates$`Ceramic Ware identification`)

#select just the dates from the Atriz valley
atriz.dates <- filter(alldates, Municipality %in% c("Pasto"))

atriz.dates <-  atriz.dates %>%
  mutate(`Ceramic Ware identification` = case_when(
    `Ceramic Ware identification` == "Piartal" ~ "1",
    `Ceramic Ware identification` == "Tuza" ~ "2",
    `Ceramic Ware identification` == "Undecorated" ~ "3",
    TRUE ~ as.character(`Ceramic Ware identification`)  # Keep all other values as they are
  ))
atriz.dates$`Ceramic Ware identification` <- as.factor(atriz.dates$`Ceramic Ware identification`)

## **Calibrate the selected dates**
atriz <- calibrate(x = atriz.dates$`Date (BP)`, errors = atriz.dates$`error (BP)`, 
                ids = atriz.dates$`radiocarbon codes`, dateDetails = atriz.dates$`Ceramic Ware identification`)
plot(atriz)

##**Plot the calibrated dates**

multiplot(atriz, calendar = "BP", label.offset = 300, decreasing = TRUE,
          gapFactor = 0.1, col.fill = atriz$metadata$Details)


pdf("Figure27.pdf", width = 11, height = 8) 
multiplot(atriz, calendar = "BP", label.offset = 300, decreasing = TRUE,
          gapFactor = 0.1, col.fill = atriz$metadata$Details)
dev.off()

# colors adjusted to follow conventions in Inkscape

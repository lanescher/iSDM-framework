

## ---------------------------
## This code was written by: C. Lane Scher
## For questions: cls7052@psu.edu
## Date Created: 2025-01-31
## ---------------------------


## ---------------------------
## Objective: Analyze output
##
## 
## Input: output from models
##   
##
## Output: Figures and summary
##
##  
## ---------------------------




library(sf)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(SpFut.flexiSDM)

na <- rnaturalearth::ne_countries(continent = "North America", 
                                  returnclass = "sf", 
                                  scale = 10) %>%
  st_transform(crs = 3857)

wa <- rnaturalearth::ne_download(type = "lakes", 
                                 category = "physical", load = T,
                                 returnclass = "sf",
                                 scale = 10) %>%
  st_transform(crs = 3857)

st <- rnaturalearth::ne_states(country = c("Canada", "Mexico", "United States of America"),
                               returnclass = "sf") %>%
  st_transform(crs = 3857)



# Figure SX: AUC ----

auc <- c()


out.dir <- "outputs/1_RACA_PNW/"
# Load data from all model runs
for (i in 1:3) {
  load(paste0(out.dir, "data", i, "-info.rdata"))
  all.auc$block <- as.character(all.auc$block)
  auc <- bind_rows(auc, all.auc)
}
# now add full model
load(paste0(out.dir, "datafull-info.rdata"))
auc <- bind_rows(auc, all.auc)

out.dir <- "outputs/2_PLSE_SE/"
# Load data from all model runs
for (i in 1:3) {
  load(paste0(out.dir, "data", i, "-info.rdata"))
  all.auc$block <- as.character(all.auc$block)
  auc <- bind_rows(auc, all.auc)
}
# now add full model
load(paste0(out.dir, "datafull-info.rdata"))
auc <- bind_rows(auc, all.auc)

out.dir <- "../hpcfiles/3_EBIS_NE/"
# Load data from all model runs
for (i in 1:3) {
  load(paste0(out.dir, "data", i, "-info.rdata"))
  all.auc$block <- as.character(all.auc$block)
  auc <- bind_rows(auc, all.auc)
}
# now add full model
load(paste0(out.dir, "datafull-info.rdata"))
auc <- bind_rows(auc, all.auc)

mean(auc$AUCout, na.rm = T)
mean(auc$AUCin, na.rm = T)

auc1 <- auc %>%
  select(sp.code, block, AUCin.full, AUCout.full) %>%
  distinct()

mean(auc1$AUCout.full, na.rm = T)
mean(auc1$AUCin.full, na.rm = T)



blockcols <- c("none" = "black", "1" = "#e79f1e", "2" = "#009e73", "3" = "#cb79a8")


auc1 <- pivot_longer(auc, cols = !c("sp.code", "block", "source")) %>%
  mutate(inout = case_when(name %in% c("AUCin", "AUCin.full", "in.cell", "in.full.cell", "in.full.n", "in.n") ~ "In sample",
                                         T ~ "Out of sample"),
                type = case_when(name %in% c("AUCin", "AUCin.full", "AUCout", "AUCout.full") ~ "AUC",
                                        name %in% c("in.cell", "in.full.cell", "out.cell", "out.full.cell") ~ "cells",
                                        name %in% c("in.full.n", "in.n", "out.full.n", "out.n") ~ "samples"),
                full = case_when(name %in% c("AUCin.full", "AUCout.full", "in.full.cell", "in.full.n", "out.full.cell", "out.full.n") ~ "full",
                                        T ~ "dataset")) %>%
  select(!name) %>%
  pivot_wider(names_from = c(type), values_from = value) %>%
  mutate(source = case_when(full == "dataset" ~ source,
                                          T ~ "All"),
         species = case_when(sp.code == "EBIS" ~ "Northern Two-lined Salamander",
                             sp.code == "PLSE" ~ "Southern Red-backed Salamander",
                             sp.code == "RACA" ~ "Cascades Frog")) %>%
  distinct()

pl <- ggplot(auc1) +
  geom_line(aes(x = source, y = AUC,
                                  group = interaction(source, block), color = as.factor(block)),
                     position = position_dodge(width = 0.6)) +
  geom_point(aes(x = source, y = AUC, shape = inout,
                                   color = as.factor(block), group = as.factor(block), size = cells),
                      position = position_dodge(width = 0.6)) +
  theme_bw() +
  scale_shape_manual(values = c(16, 1)) +
  scale_color_manual(values = blockcols) +
  labs(x = "Dataset", y = "AUC", color = "Excluded block",
                shape = "Validation", size = "Number of cells \nwith validation data") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background = element_blank(),
        legend.position = "bottom") +
  scale_x_discrete(labels = scales::label_wrap(15)) +
  facet_wrap(~species, scales = "free_x") +
  guides(color = guide_legend(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5),
         shape = guide_legend(title.position="top", title.hjust = 0.5))

ggsave(pl, file = "outputs/figures/FigSX-AUC.jpg", height = 5, width = 11)






# RACA ----
out.dir <- "outputs/1_RACA_PNW/"
load(paste0(out.dir, "datafull-info.rdata"))
load(paste0(out.dir, "datafull.rdata"))
load(paste0(out.dir, "region.rdata"))


# Load data from all model runs
auc <- c()
proc <- c()
obs <- c()
alpha <- c()
for (i in 1:3) {
  load(paste0(out.dir, "data", i, "-info.rdata"))
  load(paste0(out.dir, "data", i, ".rdata"))
  
  out$process.coef$block.out <- as.character(out$process.coef$block.out)
  out$obs.coef$block.out <- as.character(out$obs.coef$block.out)
  out$alpha$block.out <- as.character(out$alpha$block.out)
  
  auc <- bind_rows(auc, all.auc)
  proc <- bind_rows(proc, out$process.coef)
  obs <- bind_rows(obs, out$obs.coef)
  alpha <- bind_rows(alpha, out$alpha)
  
}

# now add full model
load(paste0(out.dir, "datafull-info.rdata"))

out$process.coef$block.out <- as.character(out$process.coef$block.out)
out$obs.coef$block.out <- as.character(out$obs.coef$block.out)
out$alpha$block.out <- as.character(out$alpha$block.out)

proc <- bind_rows(proc, out$process.coef)
obs <- bind_rows(obs, out$obs.coef)
alpha <- bind_rows(alpha, out$alpha)

auc$block <- as.character(auc$block)
auc <- bind_rows(auc, all.auc)



## Figure 2: range (a), intensity (b), suitable habitat (c), uncertainty (d) ----
load(paste0(out.dir, "datafull.rdata"))


bb <- st_bbox(region$region)
xlim <- c(bb[1], bb[3])
ylim <- c(bb[2], bb[4])

base <- ggplot() +
  geom_sf(data = na, fill = "gray90") +
  geom_sf(data = st, fill = "gray90", color= "gray40") +
  geom_sf(data = wa, fill = "lightsteelblue") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "lightsteelblue"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1,
                                   hjust = 1),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        title = element_text(hjust = 0.5))

mapdata <- base +
  geom_sf(data = region$region, fill = "gray85", color = NA) +
  geom_sf(data = region$range, color = "gray20", fill = NA) +
  geom_sf(data = species.data$locs$cont, aes(shape = data.type, color = source)) +
  geom_sf(data = st, fill = NA, color = "gray40") +
  geom_sf(data = wa, fill = "lightsteelblue") +
  scale_shape_manual(values = c("PO" = 4, "DND" = 0, "count" = 15, "CMR" = 17),
                     guide = "none") +
  scale_color_discrete(labels = function(x) str_wrap(x, width = 20)) +
  scale_x_continuous(breaks = c(-124, -122, -120)) +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme(legend.position = "bottom",
        legend.key = element_rect(fill = NA)) +
  guides(color = guide_legend(title.position = "top",
                              ncol = 2)) +
  labs(color = "Data source")


tmp <- region$sp.grid %>% 
  full_join(out$lambda0, by = "conus.grid.id")
q99 <- quantile(tmp$mean, 0.99, na.rm = T)
tmp$mean[which(tmp$mean > q99)] <- q99
mapint <- base +
  geom_sf(data = region$region, fill = "gray85", color = NA) +
  geom_sf(data = region$range, color = "gray20", fill = NA) +
  geom_sf(data = tmp, aes(fill = mean, color = mean)) +
  geom_sf(data = st, fill = NA, color = "gray40") +
  geom_sf(data = wa, fill = "lightsteelblue") +
  scale_x_continuous(breaks = c(-124, -122, -120)) +
  coord_sf(xlim = xlim, ylim = ylim) +
  viridis::scale_fill_viridis(guide = guide_colorbar(), 
                              option = "magma",
                              na.value = NA,
                              direction = -1,
                              breaks = c(min(tmp$mean), max(tmp$mean)),
                              labels = c("Low", "High")) +
  viridis::scale_color_viridis(guide = guide_colorbar(), 
                               option = "magma",
                               na.value = NA,
                               direction = -1,
                               breaks = c(min(tmp$mean), max(tmp$mean)),
                               labels = c("Low", "High")) +
  theme(legend.position = "bottom",
        axis.text.y = element_blank()) +
  labs(fill = "", color = "")


tmp <- region$sp.grid %>%
  full_join(out$psi0, by = "conus.grid.id") %>%
  mutate(pres = case_when(mean > 0.5 ~ "Suitable",
                          mean <= 0.5 ~ "Unsuitable"))
q99 <- quantile(tmp$mean, 0.99, na.rm = T)
tmp$mean[which(tmp$mean > q99)] <- q99
mapocc <- base +
  geom_sf(data = region$region, fill = "gray85", color = NA) +
  geom_sf(data = region$range, color = "gray20", fill = NA) +
  geom_sf(data = tmp, aes(fill = pres, color = pres)) +
  geom_sf(data = st, fill = NA, color = "gray40") +
  geom_sf(data = wa, fill = "lightsteelblue") +
  scale_x_continuous(breaks = c(-124, -122, -120)) +
  coord_sf(xlim = xlim, ylim = ylim) +
  scale_fill_manual(values = c("Unsuitable" = "white", "Suitable" = "darkblue"), na.value = NA) +
  scale_color_manual(values = c("Unsuitable" = "white", "Suitable" = "darkblue"), na.value = NA) +
  theme(legend.position = "bottom",
        axis.text.y = element_blank()) +
  labs(fill = "", color = "")

tmp <- region$sp.grid %>%
  full_join(out$lambda0, by = "conus.grid.id")
q99 <- quantile(tmp$unc.rel, 0.99, na.rm = T)
tmp$unc.rel[which(tmp$unc.rel > q99)] <- q99
mapunc <- base +
  geom_sf(data = region$region, fill = "gray85", color = NA) +
  geom_sf(data = region$range, color = "gray20", fill = NA) +
  geom_sf(data = tmp, aes(fill = unc.rel, color = unc.rel)) +
  geom_sf(data = st, fill = NA, color = "gray40") +
  geom_sf(data = wa, fill = "lightsteelblue") +
  scale_x_continuous(breaks = c(-124, -122, -120)) +
  coord_sf(xlim = xlim, ylim = ylim) +
  viridis::scale_fill_viridis(guide = guide_colorbar(), 
                              option = "magma",
                              na.value = NA,
                              direction = -1,
                              breaks = c(min(tmp$unc.rel), max(tmp$unc.rel)),
                              labels = c("Low", "High")) +
  viridis::scale_color_viridis(guide = guide_colorbar(), 
                               option = "magma",
                               na.value = NA,
                               direction = -1,
                               breaks = c(min(tmp$unc.rel), max(tmp$unc.rel)),
                               labels = c("Low", "High")) +
  theme(legend.position = "bottom",
        axis.text.y = element_blank()) +
  labs(fill = "", color = "")

maprange <- base +
  geom_sf(data = region$region, fill = "gray85", color = NA) +
  geom_sf(data = region$range, fill = NA, color = "gray35", linewidth = 0.4, aes(linetype = range.name)) +
  geom_sf(data = st, fill = NA, color = "gray40") +
  geom_sf(data = wa, fill = "lightsteelblue") +
  scale_x_continuous(breaks = c(-124, -122, -120)) +
  coord_sf(xlim = xlim, ylim = ylim) +
  viridis::scale_fill_viridis(guide = guide_colorbar(), 
                              option = "magma",
                              na.value = NA,
                              direction = -1) +
  viridis::scale_color_viridis(guide = guide_colorbar(), 
                               option = "magma",
                               na.value = NA,
                               direction = -1) +
  theme(legend.position = "bottom",
        legend.key = element_rect(fill = "white")) +
  labs(fill = "", color = "", linetype = "")

racamap <- maprange | mapint | mapocc | mapunc
racamap <- racamap + 
  plot_annotation(tag_levels = c("a"))
ggsave(racamap, file = "outputs/figures/Fig2-racamap.jpg",
       height = 7, width = 12)





## Figure 3: parameters ----
load("outputs/1_RACA_PNW/datafull-info.rdata")
cov.labs <- read.csv("data/covariate-labels.csv")

dat <- out$process.coef %>%
  mutate(x = covariate) %>%
  mutate(cov1 = gsub("2", "", covariate),
         tmp = gsub("_x_.*", "", covariate),
         quad = case_when(substr(tmp, nchar(tmp), nchar(tmp)) == 2 ~ "^2",
                          T ~ "")) %>%
  left_join(cov.labs, by = c("cov1" = "covariate")) %>%
  mutate(x = paste0(Label, quad)) %>%
  select(!tmp)

pars <- ggplot(dat) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0,
                                   hjust = 0.5,
                                   vjust = 1)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  geom_point(aes(x = x, y = mean), col='black') +
  geom_segment(aes(x = x, xend = x,
                   y = lo, yend = hi), col='black') +
  labs(y = 'Estimate', x = "Covariate")



breaks <- 0.01
# Start with process covariates
ball <- out$process.coef

# Get names of covariates, excluding intercepts, squares, and interactions
covs <- colnames(data$Xz)
rm <- grep("Z1|2|_x_", covs)
if (length(rm) > 0) covs <- covs[-rm]


# For each covariate, produce marginal effect curve
all <- c()
all.labs <- c()
for (n in 1:length(covs)) {
  cov <- covs[n]
  
  # get index of main effect
  ind <- grep(cov, ball$covariate, value = T)
  ind1 <- grep(cov, ball$covariate)
  ind1 <- ind1[ind == covs[n]]
  
  # get index of quadratic (if it exists)
  ind <- grep(paste0(cov, "2"), ball$covariate, value = T)
  ind2 <- grep(paste0(cov, "2"), ball$covariate)
  ind2 <- ind2[ind == paste0(cov, "2")]
  
  # get scaled covariate values
  q99 <- quantile(data$Xz[,cov], 0.99)
  data$Xz[,cov][data$Xz[,cov] > q99] <- q99
  s <- seq(min(data$Xz[,cov]), max(data$Xz[,cov]), breaks)
  
  
  
  # linear term with no interaction
  # b1 * x
  if (length(ind2) == 0) {
    
    b <- as.numeric(ball[ind1, "mean"])
    blo <- as.numeric(ball[ind1, "lo"])
    bhi <- as.numeric(ball[ind1, "hi"])
    
    use <- data.frame(cov = cov,
                      x = s,
                      mean = exp(b * s),
                      hi = exp(bhi * s),
                      lo = exp(blo * s),
                      factor = "none")
    
    all <- bind_rows(all, use)
  }
  
  
  # quadratic term
  # b1 * x + b2 * x^2
  if (length(ind2) == 1) {
    
    quad <- ind2
    main <- ind1
    
    b1 <- as.numeric(ball[main, "mean"])
    blo1 <- as.numeric(ball[main, "lo"])
    bhi1 <- as.numeric(ball[main, "hi"])
    
    b2 <- as.numeric(ball[quad, "mean"])
    blo2 <- as.numeric(ball[quad, "lo"])
    bhi2 <- as.numeric(ball[quad, "hi"])
    
    
    use <- data.frame(cov = cov,
                      x = s,
                      mean = exp(b1 * s + b2 * s^2),
                      hi = exp(bhi1 * s + bhi2 * s^2),
                      lo = exp(blo1 * s + blo2 * s^2),
                      factor = "none")
    
    all <- bind_rows(all, use)
  }
  
  
}


all1 <- inner_join(all, cov.labs, by = c("cov" = "covariate")) %>%
  mutate(cov = Label)


effects <- ggplot(filter(all1)) +
  geom_hline(yintercept = 0) +
  geom_ribbon(aes(x = x, ymax = hi, ymin = lo), alpha = 0.5) +
  geom_line(aes(x = x, y = mean)) +
  facet_wrap(~ cov, scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank())+
  labs(x = "Scaled covariate value", y = "Exp(Estimate)")

covs <- pars / effects
covs <- covs + 
  plot_annotation(tag_levels = "a") + 
  plot_layout(widths = c(1, 1.2))

ggsave(covs, file = "outputs/figures/Fig3-racaparameters.jpg",
       height = 6, width = 7)






# EBIS ----


load("../hpcfiles/3_EBIS_NE/datafull-info.rdata")
load("../hpcfiles/3_EBIS_NE/datafull.rdata")
load("../hpcfiles/3_EBIS_NE/region.rdata")




## Figure 4: PO data (a and b), intensity (c) ----
codeKey <- read.csv("data/model-specieslist.csv")
sp.code.all <- codeKey %>%
  filter(DS.code == "EBIS") %>%
  pull(all.codes)

allfiles <- read.csv("data/dataset-summary-full.csv")
allfiles <- allfiles[grep("EBIS", allfiles$species),] %>%
  select(-species, -percentdet) %>%
  distinct()

# detection covariates for each of these datasets
covs <- read.csv("data/00-data-summary-flexiSDM.csv") %>%
  filter(Data.Swamp.file.name %in% allfiles$file) %>%
  select(Data.Swamp.file.name, Covar.mean, Covar.sum)
covariates <- list()
for (i in 1:nrow(covs)) {
  covs.mean <- unlist(strsplit(covs$Covar.mean[i], split = ", "))
  covs.sum <- unlist(strsplit(covs$Covar.sum[i], split = ", "))
  #area <- unlist(strsplit(covs$Area[i], split = ","))
  covs1 <- c(covs.mean, covs.sum)
  covs1 <- covs1[which(is.na(covs1) == F)]
  covariates[[covs$Data.Swamp.file.name[i]]] <- covs1
}

species.data <- load_species_data("EBIS",
                                  sp.code.all,
                                  file.name = allfiles$file,
                                  file.label = allfiles$name,
                                  file.path = "data/data-ready/",
                                  keep.cols = covariates,
                                  region = region, 
                                  filter.region = T,
                                  year.start = 1994,
                                  year.end = 2025,
                                  coordunc = 1000,
                                  coordunc_na.rm = T,
                                  spat.thin = T,
                                  keep.conus.grid.id = gridkey$conus.grid.id[which(gridkey$group == "train")])


bb <- st_bbox(region$region)
xlim <- c(bb[1], bb[3])
ylim <- c(bb[2], bb[4])

base <- ggplot() +
  geom_sf(data = na, fill = "gray90") +
  geom_sf(data = st, fill = "gray90", color= "gray40") +
  geom_sf(data = wa, fill = "lightsteelblue") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "lightsteelblue"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1,
                                   hjust = 1),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.background = element_blank(),
        strip.text.x = element_blank())


use1 <- species.data$locs$cont %>%
  filter(source %in% c("VT PO", "NY Atlas", "WV PO")) %>%
  mutate(map = "state")

pl1 <- base +
  geom_sf(data = use1, aes(color = source), alpha = 0.5, size = 1) +
  geom_sf(data = region$range, color = "black", fill = NA) +
  geom_sf(data = wa, fill = "lightsteelblue") +
  coord_sf(xlim = xlim, ylim = ylim) +
  labs(color = "Source", shape = "Data type") +
  scale_color_manual(values = c("NY Atlas" = "purple",
                                "VT PO" = "blue",
                                "WV PO" = "orange")) +
  theme(legend.key = element_rect(fill = NA),
        legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3), 
                              ncol = 1,
                              title.position = "left"))

use2 <- species.data$locs$cont %>%
  filter(source %in% c("Museum", "iNaturalist")) %>%
  mutate(map = "state")

pl2 <- base +
  geom_sf(data = use2, aes(color = source), alpha = 0.5, size = 1) +
  geom_sf(data = region$range, color = "black", fill = NA) +
  geom_sf(data = wa, fill = "lightsteelblue") +
  coord_sf(xlim = xlim, ylim = ylim) +
  labs(color = "Source", shape = "Data type") +
  scale_color_manual(values = c("iNaturalist" = "cadetblue",
                                "Museum" = "darkred")) +
  theme(legend.key = element_rect(fill = NA),
        legend.position = "bottom",
        axis.text.y = element_blank()) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3), 
                              ncol = 1,
                              title.position = "left"))



lam <- out$lambda0 %>%
  full_join(region$sp.grid, ., by = "conus.grid.id")
q99 <- quantile(lam$mean, 0.99, na.rm = T)
lam$mean[which(lam$mean > q99)] <- q99

pl3 <- base +
  geom_sf(data = lam, aes(fill = mean, color = mean)) +
  #geom_sf(data = st, fill = NA, color= "gray40") +
  coord_sf(xlim = xlim, ylim = ylim) +
  labs(fill = "Relative abundance", color = "Relative abundance") +
  theme(axis.text.y = element_blank(), legend.title = element_text(hjust = 0.5)) +
  viridis::scale_fill_viridis(option = "magma",
                              na.value = "black",
                              direction = -1,
                              limits = range(c(lam$mean, lam$mean)),
                              breaks = c(min(lam$mean), max(lam$mean)),
                              labels = c("Low", "High")) +
  viridis::scale_color_viridis(option = "magma",
                               na.value = "black",
                               direction = -1,
                               limits = range(c(lam$mean, lam$mean)),
                               breaks = c(min(lam$mean), max(lam$mean)),
                               labels = c("Low", "High")) +
  guides(color = guide_colorbar(title.position = "top",
                                position = "bottom"),
         fill = guide_colorbar(title.position = "top",
                               position = "bottom"))


pl <- pl1 + pl2 + pl3

ggsave(pl, file = "outputs/figures/Fig4-ebismap.jpg",
       height = 5, width = 12)







## Figure 5: parameters ----


proc <- ggplot(filter(out$process.coef)) +
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(x = covariate, y = mean, ymin = lo, ymax = hi),
                  position = position_dodge(width = 0.3)) +
  scale_x_discrete(labels = c("elevation" = "Elevation",
                              "elevation2" = "Elevation^2",
                              "forest" = "Forest cover",
                              "prec" = "Annual \nprecipitation",
                              "streamLength.km" = "Total stream \nlength")) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 10),
        axis.text.x = element_text(vjust = 1,
                                   hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  labs(x = "Covariate", y = "Estimate")

alpha <- ggplot(filter(out$alpha)) +
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(x = name, y = mean, ymin = lo, ymax = hi),
                  position = position_dodge(width = 0.3)) +
  scale_x_discrete(labels = c("iNaturalist" = "iNaturalist",
                              "MA VES" = "MA VES",
                              "MDMBSS" = "MDMBSS",
                              "Museum" = "Museum",
                              "NY Atlas, WV PO, VT PO" = "NY Atlas, \nWV PO, \nVT PO",
                              "USGS VES - NE" = "USGS VES")) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 10),
        axis.text.x = element_text(vjust = 1,
                                   hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  labs(x = "Dataset", y = "Estimate")


state <- ggplot(filter(out$obs.coef, name %in% c("NY Atlas, WV PO, VT PO"))) +
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(x = covariate, y = mean, ymin = lo, ymax = hi),
                  position = position_dodge(width = 0.3)) +
  scale_x_discrete(labels = c("traveltime" = "Travel time \nto city",
                                "VT" = "Vermont",
                                "WV" = "West Virginia")) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 10),
        axis.text.x = element_text(vjust = 1,
                                   hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  labs(x = "Effort covariate", y = "Estimate")


pars <- proc / alpha / state
pars <- pars + 
  plot_annotation(tag_levels = c("a"))

ggsave(pars, file = "outputs/figures/Fig5-ebisparameters.jpg",
       height = 10, width = 6)




# AUC
for (i in 1:3) {
  # load model that accounts for state effort
  load(paste0("outputs/03-species-models/MVPv1/29_EBIS_proj08/data", i, "-info.rdata"))
  load("outputs/03-species-models/MVPv1/29_EBIS_proj08/datafull.rdata")
  load("outputs/03-species-models/MVPv1/29_EBIS_proj08/region.rdata")
  
}





# PLSE ----



out.dir <- "outputs/2_PLSE_SE/"
load("outputs/2_PLSE_SE/datafull-info.rdata")
load("outputs/2_PLSE_SE/datafull.rdata")
load("outputs/2_PLSE_SE/region.rdata")



# # Load data from all model runs
# auc <- c()
# proc <- c()
# obs <- c()
# alpha <- c()
# for (i in 1:3) {
#   load(paste0(out.dir, "data", i, "-info.rdata"))
#   load(paste0(out.dir, "data", i, ".rdata"))
#   
#   out$process.coef$block.out <- as.character(out$process.coef$block.out)
#   out$obs.coef$block.out <- as.character(out$obs.coef$block.out)
#   out$alpha$block.out <- as.character(out$alpha$block.out)
#   
#   auc <- bind_rows(auc, all.auc)
#   proc <- bind_rows(proc, out$process.coef)
#   obs <- bind_rows(obs, out$obs.coef)
#   alpha <- bind_rows(alpha, out$alpha)
#   
# }
# 
# # now add full model
# load(paste0(out.dir, "datafull-info.rdata"))
# 
# out$process.coef$block.out <- as.character(out$process.coef$block.out)
# out$obs.coef$block.out <- as.character(out$obs.coef$block.out)
# out$alpha$block.out <- as.character(out$alpha$block.out)
# 
# proc <- bind_rows(proc, out$process.coef)
# obs <- bind_rows(obs, out$obs.coef)
# alpha <- bind_rows(alpha, out$alpha)
# 
# 
# auc$block <- as.character(auc$block)
# auc <- bind_rows(auc, all.auc)

## Figure 6: map ----

bb <- st_bbox(region$region)
xlim <- c(bb[1], bb[3])
ylim <- c(bb[2], bb[4])


base <- ggplot() +
  geom_sf(data = na, fill = "gray90") +
  geom_sf(data = st, fill = "gray90", color= "gray40") +
  geom_sf(data = wa, fill = "lightsteelblue") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "lightsteelblue"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1,
                                   hjust = 1),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        title = element_text(hjust = 0.5))




maprange <- base +
  geom_sf(data = region$region, fill = "gray85", color = NA) +
  geom_sf(data = region$range, fill = NA, color = "gray35", linewidth = 0.4, aes(linetype = range.name)) +
  geom_sf(data = st, fill = NA, color = "gray40") +
  geom_sf(data = wa, fill = "lightsteelblue") +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme(legend.position = "top",
        legend.key = element_rect(fill = NA),
        axis.text.x = element_blank()) +
  guides(color = guide_legend(title.position = "top",
                              ncol = 2)) +
  labs(color = "Data source", linetype = "")



tmp <- out$lambda0 %>%
  full_join(region$sp.grid, ., by = "conus.grid.id")
q99 <- quantile(tmp$mean, 0.99, na.rm = T)
tmp$mean[which(tmp$mean > q99)] <- q99


mapint <- base +
  geom_sf(data = region$region, fill = "gray85", color = NA) +
  geom_sf(data = region$range, color = "gray20", fill = NA) +
  geom_sf(data = tmp, aes(fill = mean, color = mean)) +
  geom_sf(data = st, fill = NA, color = "gray40") +
  geom_sf(data = wa, fill = "lightsteelblue") +
  coord_sf(xlim = xlim, ylim = ylim) +
  viridis::scale_fill_viridis(guide = guide_colorbar(), 
                              option = "magma",
                              na.value = NA,
                              direction = -1,
                              breaks = c(min(tmp$mean), max(tmp$mean)),
                              labels = c("Low", "High")) +
  viridis::scale_color_viridis(guide = guide_colorbar(), 
                               option = "magma",
                               na.value = NA,
                               direction = -1,
                               breaks = c(min(tmp$mean), max(tmp$mean)),
                               labels = c("Low", "High")) +
  theme(legend.position = "top",
        axis.text.y = element_blank(),
        axis.text.x = element_blank()) +
  labs(fill = "", color = "")


tmp <- out$XB0 %>%
  full_join(region$sp.grid, ., by = "conus.grid.id")
q99 <- quantile(tmp$mean, 0.99, na.rm = T)
tmp$mean[which(tmp$mean > q99)] <- q99
tmp$mean <- exp(tmp$mean)


mapxb <- base +
  geom_sf(data = region$region, fill = "gray85", color = NA) +
  geom_sf(data = region$range, color = "gray20", fill = NA) +
  geom_sf(data = tmp, aes(fill = mean, color = mean)) +
  geom_sf(data = st, fill = NA, color = "gray40") +
  geom_sf(data = wa, fill = "lightsteelblue") +
  coord_sf(xlim = xlim, ylim = ylim) +
  viridis::scale_fill_viridis(guide = guide_colorbar(), 
                              option = "magma",
                              na.value = NA,
                              direction = -1,
                              breaks = c(min(tmp$mean), max(tmp$mean)),
                              labels = c("Low", "High")) +
  viridis::scale_color_viridis(guide = guide_colorbar(), 
                               option = "magma",
                               na.value = NA,
                               direction = -1,
                               breaks = c(min(tmp$mean), max(tmp$mean)),
                               labels = c("Low", "High")) +
  theme(legend.position = "bottom") +
  labs(fill = "", color = "")


tmp <- out$spat %>%
  full_join(region$sp.grid, ., by = "conus.grid.id")


mapspat <- base +
  geom_sf(data = region$region, fill = "gray85", color = NA) +
  geom_sf(data = region$range, color = "gray20", fill = NA) +
  geom_sf(data = tmp, aes(fill = mean, color = mean)) +
  geom_sf(data = st, fill = NA, color = "gray40") +
  geom_sf(data = wa, fill = "lightsteelblue") +
  coord_sf(xlim = xlim, ylim = ylim) +
  scale_fill_gradient2(guide = guide_colorbar(), high = "darkblue", low = "darkred", mid = "white", na.value = NA) +
  scale_color_gradient2(guide = "none", high = "darkblue", low = "darkred", mid = "white", na.value = NA) +
  theme(legend.position = "bottom",
        axis.text.y = element_blank()) +
  labs(fill = "", color = "")



plsemap <- (maprange | mapint) / (mapxb | mapspat)
plsemap <- plsemap + 
  plot_annotation(tag_levels = c("a"))
ggsave(plsemap, file = "outputs/figures/Fig6-plsemap.jpg",
       height = 7, width = 7)


## Figure 7: parameters ----
load("outputs/2_PLSE_SE/datafull-info.rdata")
cov.labs <- read.csv("data/covariate-labels.csv")

dat <- out$process.coef %>%
  mutate(x = covariate) %>%
  mutate(cov1 = gsub("2", "", covariate),
         tmp = gsub("_x_.*", "", covariate),
         quad = case_when(substr(tmp, nchar(tmp), nchar(tmp)) == 2 ~ "^2",
                          T ~ "")) %>%
  left_join(cov.labs, by = c("cov1" = "covariate")) %>%
  mutate(x = paste0(Label, quad)) %>%
  select(!tmp)

pars <- ggplot(dat) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0,
                                   hjust = 0.5,
                                   vjust = 1)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  geom_point(aes(x = x, y = mean), col='black') +
  geom_segment(aes(x = x, xend = x,
                   y = lo, yend = hi), col='black') +
  labs(y = 'Estimate', x = "Covariate")



breaks <- 0.01
# Start with process covariates
ball <- out$process.coef

# Get names of covariates, excluding intercepts, squares, and interactions
covs <- colnames(data$Xz)
rm <- grep("Z1|2|_x_", covs)
if (length(rm) > 0) covs <- covs[-rm]


# For each covariate, produce marginal effect curve
all <- c()
all.labs <- c()
for (n in 1:length(covs)) {
  cov <- covs[n]
  
  # get index of main effect
  ind <- grep(cov, ball$covariate, value = T)
  ind1 <- grep(cov, ball$covariate)
  ind1 <- ind1[ind == covs[n]]
  
  # get index of quadratic (if it exists)
  ind <- grep(paste0(cov, "2"), ball$covariate, value = T)
  ind2 <- grep(paste0(cov, "2"), ball$covariate)
  ind2 <- ind2[ind == paste0(cov, "2")]
  
  # get scaled covariate values
  q99 <- quantile(data$Xz[,cov], 0.99)
  data$Xz[,cov][data$Xz[,cov] > q99] <- q99
  s <- seq(min(data$Xz[,cov]), max(data$Xz[,cov]), breaks)
  
  
  
  # linear term with no interaction
  # b1 * x
  if (length(ind2) == 0) {
    
    b <- as.numeric(ball[ind1, "mean"])
    blo <- as.numeric(ball[ind1, "lo"])
    bhi <- as.numeric(ball[ind1, "hi"])
    
    use <- data.frame(cov = cov,
                      x = s,
                      mean = exp(b * s),
                      hi = exp(bhi * s),
                      lo = exp(blo * s),
                      factor = "none")
    
    all <- bind_rows(all, use)
  }
  
  
  # quadratic term
  # b1 * x + b2 * x^2
  if (length(ind2) == 1) {
    
    quad <- ind2
    main <- ind1
    
    b1 <- as.numeric(ball[main, "mean"])
    blo1 <- as.numeric(ball[main, "lo"])
    bhi1 <- as.numeric(ball[main, "hi"])
    
    b2 <- as.numeric(ball[quad, "mean"])
    blo2 <- as.numeric(ball[quad, "lo"])
    bhi2 <- as.numeric(ball[quad, "hi"])
    
    
    use <- data.frame(cov = cov,
                      x = s,
                      mean = exp(b1 * s + b2 * s^2),
                      hi = exp(bhi1 * s + bhi2 * s^2),
                      lo = exp(blo1 * s + blo2 * s^2),
                      factor = "none")
    
    all <- bind_rows(all, use)
  }
  
  
}


all1 <- inner_join(all, cov.labs, by = c("cov" = "covariate")) %>%
  mutate(cov = Label)


effects <- ggplot(filter(all1)) +
  geom_hline(yintercept = 0) +
  geom_ribbon(aes(x = x, ymax = hi, ymin = lo), alpha = 0.5) +
  geom_line(aes(x = x, y = mean)) +
  facet_wrap(~ cov, scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank())+
  labs(x = "Scaled covariate value", y = "Exp(Estimate)")

covs <- pars / effects
covs <- covs + 
  plot_annotation(tag_levels = "a") + 
  plot_layout(widths = c(1, 1.2))

ggsave(covs, file = "outputs/figures/Fig7-plseparameters.jpg",
       height = 6, width = 7)





# Table SX: Data summary ----

dat <- c()

load("outputs/1_RACA_PNW/datafull-info.rdata")
dat1 <- species.data$locs$disc %>%
  mutate(species = "Cascades Frog",
         sp.code = "RACA")
dat <- bind_rows(dat, dat1)

load("outputs/2_PLSE_SE/datafull-info.rdata")
dat1 <- species.data$locs$disc %>%
  mutate(species = "Southern Red-backed Salamander",
         sp.code = "PLSE")
dat <- bind_rows(dat, dat1)

load("outputs/3_EBIS_NE/datafull-info.rdata")
dat1 <- species.data$locs$disc %>%
  mutate(species = "Northern Two-lined Salamander",
         sp.code = "EBIS")
dat <- bind_rows(dat, dat1)


n.obs <- dat %>%
  select(species, sp.code, source, data.type, survey.id) %>%
  distinct() %>%
  group_by(species, sp.code, source, data.type) %>%
  summarize(n.obs = n(), .groups = "drop")

n.site <- dat %>%
  select(species, sp.code, source, data.type, site.id) %>%
  distinct() %>%
  group_by(species, sp.code, source, data.type) %>%
  summarize(n.site = n(), .groups = "drop")

n.cell <- dat %>%
  select(species, sp.code, source, data.type, conus.grid.id) %>%
  distinct() %>%
  group_by(species, sp.code, source, data.type) %>%
  summarize(n.cell = n(), .groups = "drop")

datsum <- full_join(n.obs, n.site, by = c("species", "source", "data.type")) %>%
  full_join(n.cell, by = c("species", "source", "data.type")) %>%
  select(species, sp.code, data.type, source, n.obs, n.site, n.cell)


sum <- read.csv("data/00-data-summary-flexiSDM.csv")

datsum <- inner_join(datsum, sum, by = c("sp.code" = "Species",
                                         "source" = "Name")) %>%
  select(species, data.type, source, Method, n.obs, n.site, n.cell, 
         Covar.mean, Covar.sum, PO.extent) %>%
  mutate(covariates = paste0(Covar.mean, ", ", Covar.sum)) %>%
  select(-Covar.sum, -Covar.mean) %>%
  arrange(species, data.type, desc(n.obs))
datsum$covariates[substr(datsum$covariates, 1, 1) == ","] <- gsub(", ", "", datsum$covariates[substr(datsum$covariates, 1, 1) == ","])
datsum$covariates[grep("^yday, $", datsum$covariates)] <- 'yday'


colnames(datsum) <- c("Species", "Data type", "Dataset", "Method", "Number of observations", 
                      "Number of sites", "Number of hexbins", "Extent", "Covariates")

write.csv(datsum, file = "outputs/tables/TabSX-datasummary.csv",
          row.names = F)



# Get extent of each dataset

load("outputs/1_RACA_PNW/datafull-info.rdata")
load("outputs/1_RACA_PNW/region.rdata")

ggplot() + geom_sf(data = region$range) +
  geom_sf(data = species.data$locs$cont) +
  facet_wrap(~ source)

load("outputs/2_PLSE_SE/datafull-info.rdata")
load("outputs/2_PLSE_SE/region.rdata")

ggplot() + geom_sf(data = region$range) +
  geom_sf(data = species.data$locs$cont) +
  facet_wrap(~ source)


load("outputs/3_EBIS_NE/datafull-info.rdata")
load("outputs/3_EBIS_NE/region.rdata")

range <- region$range %>%
  summarize(geometry = st_union(geometry))

ggplot() + geom_sf(data = range) +
  geom_sf(data = species.data$locs$cont) +
  facet_wrap(~ source)




# Figure SX: CV maps ----
blockcols <- c("none" = "black", "1" = "#e79f1e", "2" = "#009e73", "3" = "#cb79a8")



base <- ggplot() +
  geom_sf(data = na, fill = "gray90") +
  geom_sf(data = st, fill = "gray90", color= "gray40") +
  geom_sf(data = wa, fill = "lightsteelblue") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "lightsteelblue"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1,
                                   hjust = 1),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        title = element_text(hjust = 0.5))




load(paste0("outputs/1_RACA_PNW/datafull-info.rdata"))
load("outputs/1_RACA_PNW/setup_none.Rdata")

blocks <- sb$blocks %>%
  st_intersection(region$region)

bb <- st_bbox(region$region)
xlim <- c(bb[1], bb[3])
ylim <- c(bb[2], bb[4])


raca <- base +
  #geom_sf(data = species.data$locs$cont, aes(shape = data.type, color = source)) +
  geom_sf(data = wa, fill = "lightsteelblue") +
  geom_sf(data = blocks, aes(fill = as.character(folds)), alpha = 0.5, color = NA) +
  geom_sf(data = st, fill = NA, color = "gray40") +
  #geom_sf(data = region$range, color = "gray35", fill = NA, linewidth = 0.4, aes(linetype = range.name)) +
  scale_shape_manual(values = c("PO" = 4, "DND" = 0, "count" = 15, "CMR" = 17),
                     guide = "none") +
  scale_color_discrete(labels = function(x) str_wrap(x, width = 20), guide = "none") +
  scale_fill_manual(values = blockcols, guide = "none") +
  scale_linetype(guide = "none") +
  scale_x_continuous(breaks = c(-124, -122, -120)) +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme(legend.position = "bottom",
        legend.key = element_rect(fill = NA)) +
  # guides(color = guide_legend(title.position = "top",
  #                             ncol = 2)) +
  labs(color = "Data source")



load(paste0("outputs/2_PLSE_SE/datafull-info.rdata"))
load("outputs/2_PLSE_SE/region.rdata")
#load("outputs/2_PLSE_SE/setup_none.Rdata")

# blocks <- sb$blocks %>%
#   st_intersection(region$region)


bb <- st_bbox(region$region)
xlim <- c(bb[1], bb[3])
ylim <- c(bb[2], bb[4])


plse <- base +
  #geom_sf(data = species.data$locs$cont, aes(shape = data.type, color = source)) +
  geom_sf(data = wa, fill = "lightsteelblue") +
  #geom_sf(data = blocks, aes(fill = as.character(folds)), alpha = 0.5, color = NA) +
  geom_sf(data = st, fill = NA, color = "gray40") +
  #geom_sf(data = region$range, color = "gray35", linewidth = 0.4, fill = NA, aes(linetype = range.name)) +
  scale_shape_manual(values = c("PO" = 4, "DND" = 0, "count" = 15, "CMR" = 17),
                     guide = "none") +
  scale_color_discrete(labels = function(x) str_wrap(x, width = 20), guide = "none") +
  scale_fill_manual(values = blockcols, guide = "none") +
  scale_linetype(guide = "none") +
  scale_x_continuous(breaks = c(-124, -122, -120)) +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme(legend.position = "bottom",
        legend.key = element_rect(fill = NA)) +
  # guides(color = guide_legend(title.position = "top",
  #                             ncol = 2)) +
  labs(color = "Data source")




load(paste0("outputs/3_EBIS_NE/datafull-info.rdata"))
load("outputs/3_EBIS_NE/region.rdata")
load("outputs/3_EBIS_NE/setup_none.Rdata")

bb <- st_bbox(region$region)
xlim <- c(bb[1], bb[3])
ylim <- c(bb[2], bb[4])


blocks <- sb$blocks %>%
  st_intersection(region$region)
# ggplot(blocks) + geom_sf()

ebis <- base +
  #geom_sf(data = species.data$locs$cont, aes(shape = data.type, color = source)) +
  geom_sf(data = wa, fill = "lightsteelblue") +
  geom_sf(data = blocks, aes(fill = as.character(folds)), alpha = 0.5, color = NA) +
  geom_sf(data = st, fill = NA, color = "gray40") +
  #geom_sf(data = region$range, color = "gray35", linewidth = 0.4, fill = NA, aes(linetype = range.name)) +
  scale_shape_manual(values = c("PO" = 4, "DND" = 0, "count" = 15, "CMR" = 17),
                     guide = "none") +
  scale_color_discrete(labels = function(x) str_wrap(x, width = 20), guide = "none") +
  scale_fill_manual(values = blockcols, guide = "none") +
  scale_linetype_manual(guide = "none", values = c("IUCN" = "dotted")) +
  scale_x_continuous(breaks = c(-124, -122, -120)) +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme(legend.position = "bottom",
        legend.key = element_rect(fill = NA)) +
  # guides(color = guide_legend(title.position = "top",
  #                             ncol = 2)) +
  labs(color = "Data source")



pl <- raca | ebis / plse
pl <- pl + plot_annotation(tag_levels = "a")

ggsave(pl, file = "outputs/figures/FigSX-CV.jpg", height = 9, width = 12)




# Appendix 3: tau ----

dirs <- list.dirs("outputs/", recursive = F)
dirs <- dirs[grep("RACA", dirs)]

tau <- c()
for (d in 1:length(dirs)) {
  load(paste0(dirs[d], "/datafull.rdata"))
  tau1 <- out$tau %>%
    mutate(model = d)
  
  tau <- bind_rows(tau, tau1)
}

priors <- c()
for (p in 4:6) {
  if (p == 4) vec <- rgamma(100000, 5, 5)
  if (p == 5) vec <- rnorm(100000, 0, (1/sqrt(0.1)))
  if (p == 6) vec <- rgamma(100000, 0.01, 0.01)
  
  tmp <- data.frame(model = p,
                    prior = vec)
  

  priors <- bind_rows(priors, tmp)
}

priors <- priors %>%
  mutate(name = case_when(model == 4 ~ "Gamma(5, 5)",
                          model == 5 ~ "Normal(0, 3.16)",
                          model == 6 ~ "Gamma(0.01, 0.01)"))

tau <- tau %>%
  mutate(name = case_when(model == 4 ~ "Gamma(5, 5)",
                          model == 5 ~ "Normal(0, 3.16)",
                          model == 6 ~ "Gamma(0.01, 0.01)"))

ggplot() + 
  geom_violin(data = priors, aes(x = name, y = prior)) +
  geom_pointrange(data = filter(tau, lo != hi), aes(x = name, y = mean, ymin = lo, ymax = hi)) +
  theme_bw() +
  labs(x = "Model", y = "Value") +
  facet_wrap(~model, scales = "free_x") +
  coord_cartesian(ylim = c(-10, 10))

ggplot() + 
  geom_violin(data = priors, aes(x = name, y = prior)) +
  geom_pointrange(data = filter(tau, lo != hi), aes(x = name, y = mean, ymin = lo, ymax = hi)) +
  theme_bw() +
  labs(x = "Model", y = "Value") +
  facet_wrap(~model, scales = "free_x") +
  coord_cartesian(ylim = c(-0.25, 0.25))


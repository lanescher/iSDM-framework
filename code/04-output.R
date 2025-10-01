

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



# Figure 1: Model framework ----



# out.dir <- "outputs/2_RACA_tau1//"
# load(paste0(out.dir, "datafull-info.rdata"))
# load(paste0(out.dir, "datafull.rdata"))
# load(paste0(out.dir, "region.rdata"))


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
# auc$block <- as.character(auc$block)
# auc <- bind_rows(auc, all.auc)
# 


## Figure 2: RACA- range (a), intensity (b), suitable habitat (c), uncertainty (d) ----
# load(paste0(out.dir, "datafull.rdata"))
out.dir <- "outputs/2_RACA_tau1//"
load(paste0(out.dir, "datafull-info.rdata"))
load(paste0(out.dir, "datafull.rdata"))
load(paste0(out.dir, "region.rdata"))


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
  mutate(pres = case_when(mean > 0.5 ~ "Occupied",
                          mean <= 0.5 ~ "Unoccupied"))
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
  scale_fill_manual(values = c("Unoccupied" = "white", "Occupied" = "darkblue"), na.value = NA) +
  scale_color_manual(values = c("Unoccupied" = "white", "Occupied" = "darkblue"), na.value = NA) +
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





## Figure 3: RACA- parameters ----
load("outputs/2_RACA_tau1/datafull-info.rdata")
load("outputs/2_RACA_tau1/datafull.rdata")
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
       height = 6, width = 8)






## Figure 4: GPOR- PO data and effort ----

load("outputs/8_GPOR_tau1/datafull-info.rdata")
load("outputs/8_GPOR_tau1/datafull.rdata")
load("outputs/8_GPOR_tau1/region.rdata")




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


codeKey <- read.csv("data/model-specieslist.csv")
sp.code.all <- codeKey %>%
  filter(DS.code == "GPOR") %>%
  pull(all.codes)

allfiles <- read.csv("data/dataset-summary-full.csv")
allfiles <- allfiles[grep("GPOR", allfiles$species),] %>%
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

species.data <- load_species_data("GPOR",
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
  filter(source %in% c("VT PO", "NY Atlas", "WV PO", "MA PO")) %>%
  mutate(map = "state")
use1$source <- factor(use1$source, levels = c("VT PO", "NY Atlas", "WV PO", "MA PO", "iNaturalist"))

pl1 <- base +
  geom_sf(data = use1, aes(color = source), alpha = 0.5, size = 0.5, show.legend = T) +
  geom_sf(data = region$range, color = "black", fill = NA) +
  geom_sf(data = wa, fill = "lightsteelblue") +
  coord_sf(xlim = xlim, ylim = ylim) +
  labs(color = "Source", shape = "Data type") +
  scale_color_manual(values = c("NY Atlas" = "#CC79A7",
                                "VT PO" = "#0072B2",
                                "WV PO" = "#E69F00",
                                "MA PO" = "#009E73",
                                "iNaturalist" = "#D55E00"),
                     drop = F) +
  theme(legend.key = element_rect(fill = NA),
        legend.position = "bottom",
        axis.text.x = element_blank()) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3), 
                              ncol = 1,
                              title.position = "top"))

use2 <- species.data$locs$cont %>%
  filter(source %in% c("iNaturalist")) %>%
  mutate(map = "state")
use2$source <- factor(use2$source, levels = c("VT PO", "NY Atlas", "WV PO", "MA PO", "iNaturalist"))

pl2 <- base +
  geom_sf(data = use2, aes(color = source), alpha = 0.5, size = 0.5, show.legend = T) +
  geom_sf(data = region$range, color = "black", fill = NA) +
  geom_sf(data = wa, fill = "lightsteelblue") +
  coord_sf(xlim = xlim, ylim = ylim) +
  labs(color = "Source", shape = "Data type") +
  scale_color_manual(values = c("NY Atlas" = "#CC79A7",
                                "VT PO" = "#0072B2",
                                "WV PO" = "#E69F00",
                                "MA PO" = "#009E73",
                                "iNaturalist" = "#D55E00"),
                     drop = F) +
  theme(legend.key = element_rect(fill = NA),
        legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.text.x = element_blank()) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3), 
                              ncol = 1,
                              title.position = "top"))




tmp <- out$effort %>%
  filter(PO.dataset.name %in% c("NY Atlas, WV PO, VT PO, MA PO")) %>%
  full_join(region$sp.grid, ., by = "conus.grid.id") %>%
  filter(mean != 0)

pl3 <- base +
  geom_sf(data = tmp, aes(fill = log(mean), color = log(mean))) +
  geom_sf(data = region$range, color = "black", fill = NA) +
  coord_sf(xlim = xlim, ylim = ylim) +
  viridis::scale_fill_viridis(option = "magma",
                              na.value = "black",
                              direction = -1,
                              limits = range(c(log(tmp$mean), log(tmp$mean))),
                              breaks = c(log(min(tmp$mean)), log(max(tmp$mean))),
                              labels = c("Low", "High")) +
  viridis::scale_color_viridis(option = "magma",
                               na.value = "black",
                               direction = -1,
                               limits = range(c(log(tmp$mean), log(tmp$mean))),
                               breaks = c(log(min(tmp$mean)), log(max(tmp$mean))),
                               labels = c("Low", "High")) +
  guides(fill = guide_colorbar(title.position = "top",
                               position = "right"),
         color = guide_colorbar(title.position = "top",
                                position = "right")) +
  labs(color = "Log(Effort)", fill = "Log(Effort)")


tmp <- out$effort %>%
  filter(PO.dataset.name %in% c("iNaturalist")) %>%
  full_join(region$sp.grid, ., by = "conus.grid.id") %>%
  filter(mean != 0)

pl4 <- base +
  geom_sf(data = tmp, aes(fill = log(mean), color = log(mean))) +
  geom_sf(data = region$range, color = "black", fill = NA) +
  coord_sf(xlim = xlim, ylim = ylim) +
  viridis::scale_fill_viridis(option = "magma",
                              na.value = "black",
                              direction = -1,
                              limits = range(c(log(tmp$mean), log(tmp$mean))),
                              breaks = c(log(min(tmp$mean)), log(max(tmp$mean))),
                              labels = c("Low", "High")) +
  viridis::scale_color_viridis(option = "magma",
                               na.value = "black",
                               direction = -1,
                               limits = range(c(log(tmp$mean), log(tmp$mean))),
                               breaks = c(log(min(tmp$mean)), log(max(tmp$mean))),
                               labels = c("Low", "High")) +
  guides(fill = guide_colorbar(title.position = "top",
                               position = "right"),
         color = guide_colorbar(title.position = "top",
                                position = "right")) +
  labs(color = "Log(Effort)", fill = "Log(Effort)") +
  theme(axis.text.y = element_blank())




pars1 <- out$obs.coef %>%
  filter(name %in% c("NY Atlas, WV PO, VT PO, MA PO")) %>%
  mutate(covariate = case_when(covariate == "traveltime" ~ "Travel time",
                               T ~ covariate))
pars1$covariate <- factor(pars1$covariate, levels = c("WV", "MA", "VT", "Travel time"))

pl5 <- ggplot(pars1) +
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(x = covariate, y = mean,
                      ymin = lo, ymax = hi), size = 0.2) +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  labs(x = "Covariate", y = "Estimate")

pars2 <- out$obs.coef %>%
  filter(name %in% c("iNaturalist")) %>%
  mutate(covariate = case_when(covariate == "n.inat" ~ "Number of iNaturalist records",
                               covariate == "intercept" ~ "Intercept"))

pl6 <- ggplot(pars2) +
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(x = covariate, y = mean,
                      ymin = lo, ymax = hi), size = 0.2) +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  labs(x = "Covariate", y = "")

pla <- pl1 + pl2 +
  plot_layout(guides = "collect")

plb <- pl5 + pl6

plc <- pl3 + pl4 +
  plot_layout(guides = "collect")

layout <- "
AAAA
AAAA
AAAA
BBBB
BBBB
BBBB
CCCC
"

pl <- pla / plc / plb
pl <- pl + 
  plot_annotation(tag_levels = "a") +
  plot_layout(design = layout)

ggsave(pl, file = "outputs/figures/Fig4-gporPO.jpg",
       height = 12, width = 12)







# Figure 5: GPOR- Process model ----

load("outputs/8_GPOR_tau1/datafull-info.rdata")
load("outputs/8_GPOR_tau1/datafull.rdata")
load("outputs/8_GPOR_tau1/region.rdata")

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
        axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1,
                                   hjust = 1),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 12),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        title = element_text(hjust = 0.5))

betas <- out$process.coef %>%
  mutate(covariate = case_when(covariate == "elevation" ~ "Elevation",
                               covariate == "elevation2" ~ "Elevation^2",
                               covariate == "forest" ~ "Forest",
                               covariate == "prec" ~ "Precipitation",
                               covariate == "streamLength.km" ~ "Stream length"))

plotbetas <- ggplot(betas) +
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(x = covariate, y = mean,
                      ymin = lo, ymax = hi), size = 0.2) +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  labs(x = "Covariate", y = "Estimate")




lam <- out$lambda0 %>%
  full_join(region$sp.grid, ., by = "conus.grid.id")
q99 <- quantile(lam$mean, 0.99, na.rm = T)
lam$mean[which(lam$mean > q99)] <- q99

maplam <- base +
  geom_sf(data = lam, aes(fill = mean, color = mean)) +
  #geom_sf(data = st, fill = NA, color= "gray40") +
  coord_sf(xlim = xlim, ylim = ylim) +
  labs(fill = "Relative abundance", color = "Relative abundance") +
  theme(legend.title = element_text(hjust = 0.5),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(1.5, "cm")) +
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


xb <- out$XB0 %>%
  full_join(region$sp.grid, ., by = "conus.grid.id")
q99 <- quantile(xb$mean, 0.99, na.rm = T)
xb$mean[which(xb$mean > q99)] <- q99

mapxb <- base +
  geom_sf(data = xb, aes(fill = mean, color = mean)) +
  #geom_sf(data = st, fill = NA, color= "gray40") +
  coord_sf(xlim = xlim, ylim = ylim) +
  labs(fill = "XB", color = "XB") +
  theme(axis.text.y = element_blank(), legend.title = element_text(hjust = 0.5),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(1.5, "cm")) +
  viridis::scale_fill_viridis(option = "magma",
                              na.value = "black",
                              direction = -1,
                              limits = range(c(xb$mean, xb$mean)),
                              breaks = c(min(xb$mean), max(xb$mean)),
                              labels = c("Low", "High")) +
  viridis::scale_color_viridis(option = "magma",
                               na.value = "black",
                               direction = -1,
                               limits = range(c(xb$mean, xb$mean)),
                               breaks = c(min(xb$mean), max(xb$mean)),
                               labels = c("Low", "High")) +
  guides(color = guide_colorbar(title.position = "top",
                                position = "bottom"),
         fill = guide_colorbar(title.position = "top",
                               position = "bottom"))



spat <- out$spat %>%
  full_join(region$sp.grid, ., by = "conus.grid.id")
q99 <- quantile(spat$mean, 0.99, na.rm = T)
spat$mean[which(spat$mean > q99)] <- q99

mapspat <- base +
  geom_sf(data = spat, aes(fill = mean, color = mean)) +
  #geom_sf(data = st, fill = NA, color= "gray40") +
  coord_sf(xlim = xlim, ylim = ylim) +
  labs(fill = "Spatial effect", color = "Spatial effect") +
  theme(axis.text.y = element_blank(), legend.title = element_text(hjust = 0.5),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(1.5, "cm")) +
  ggplot2::scale_fill_gradient2(guide = ggplot2::guide_colorbar(), high = "darkblue", low = "darkred", mid = "white", na.value = NA) +
  ggplot2::scale_color_gradient2(guide = "none", high = "darkblue", low = "darkred", mid = "white", na.value = NA) +
  guides(color = guide_colorbar(title.position = "top",
                                position = "bottom"),
         fill = guide_colorbar(title.position = "top",
                               position = "bottom"))




layout <- "
AAA
AAA
BBB
BBB
BBB
"
pl <- plotbetas / (maplam | mapxb | mapspat) +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = "a")

ggsave(pl, file = "outputs/figures/Fig4-gporprocess.jpg",
       height = 8, width = 12)







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




# Figure S2: CV maps ----
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
        title = element_text(hjust = 0.5))



# RACA
load("outputs/2_RACA_tau1/datafull-info.rdata")
load("outputs/2_RACA_tau1/region.rdata")

spatblocks <- make_CV_blocks(region, rows = 5, cols = 5, k = 3)


blocks <- spatblocks %>%
  st_intersection(region$region)

bb <- st_bbox(region$region)
xlim <- c(bb[1], bb[3])
ylim <- c(bb[2], bb[4])



raca <- base +
  geom_sf(data = wa, fill = "lightsteelblue") +
  geom_sf(data = blocks, aes(fill = as.character(folds)), alpha = 0.5, color = NA) +
  geom_sf(data = st, fill = NA, color = "gray40") +
  scale_fill_manual(values = blockcols) +
  scale_x_continuous(breaks = c(-124, -122, -120)) +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme(legend.key = element_rect(fill = NA)) +
  labs(fill = "Excluded\nfold")




# GPOR
load("outputs/8_GPOR_tau1/datafull-info.rdata")
load("outputs/8_GPOR_tau1/region.rdata")

spatblocks <- make_CV_blocks(region, rows = 5, cols = 5, k = 3)


blocks <- spatblocks %>%
  st_intersection(region$region)

bb <- st_bbox(region$region)
xlim <- c(bb[1], bb[3])
ylim <- c(bb[2], bb[4])


gpor <- base +
  geom_sf(data = wa, fill = "lightsteelblue") +
  geom_sf(data = blocks, aes(fill = as.character(folds)), alpha = 0.5, color = NA) +
  geom_sf(data = st, fill = NA, color = "gray40") +
  scale_fill_manual(values = blockcols) +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme(legend.key = element_rect(fill = NA)) +
  labs(fill = "Excluded\nfold")






pl <- raca | gpor
pl <- pl + 
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = 'collect')

ggsave(pl, file = "outputs/figures/FigS2-CV.jpg", height = 9, width = 12)


# Figure S3: RACA data map from auto-generated output ----



# Figure S4: GPOR data map from auto-generated output ----


# Figure S5: coarse spatial effect ----
load("outputs/8_GPOR_tau1/region.rdata")
spatRegion <- suppressWarnings(make_spatkey(region$sp.grid))


coarse <- spatRegion$spat.grid %>%
  st_transform(crs = 4326)
fine <- region$sp.grid %>%
  st_transform(crs = 4326)


main <- ggplot() +
  geom_sf(data = coarse, fill = "white", linewidth = 0.1) +
  geom_rect(aes(xmin = -86.5, xmax = -84.6, ymin = 36.6, ymax = 37.75), color = "black", fill = NA) +
  theme_bw()
inset <- ggplot() +
  geom_sf(data = fine, fill = "white", color = "gray80", linewidth = 0.4) +
  geom_sf(data = coarse, fill = NA, color = "gray20", linewidth = 0.5) +
  theme_bw() +
  coord_sf(ylim = c(36.6, 37.8), xlim = c(-86.5, -84.6), expand = F) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())


pl <- main | inset
pl <- pl + plot_annotation(tag_levels = "a")


ggsave(pl, file = "outputs/figures/FigS5-coarsegrid.jpg", height = 6, width = 12)



# Figure S6: AUC ----

blockcols <- c("none" = "black", "1" = "#e79f1e", "2" = "#009e73", "3" = "#cb79a8")

auc <- c()


out.dir <- "outputs/2_RACA_tau1/"
# Load data from all model runs
for (i in 1:3) {
  load(paste0(out.dir, "data", i, "-info.rdata"))
  all.auc$block <- as.character(all.auc$block)
  auc <- bind_rows(auc, all.auc)
}
# now add full model
load(paste0(out.dir, "datafull-info.rdata"))
auc <- bind_rows(auc, all.auc)

out.dir <- "outputs/8_GPOR_tau1/"
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
  distinct() %>%
  mutate(species = case_when(sp.code == "GPOR" ~ "Spring Salamander",
                             sp.code == "RACA" ~ "Cascades Frog"))

mean(auc1$AUCout.full, na.rm = T)
mean(auc1$AUCin.full, na.rm = T)

auc2 <- pivot_longer(auc1, cols = !c("species", "sp.code", "block"))
pl <- ggplot(auc2) +
  geom_line(aes(x = block, y = value,
                group = interaction(block), color = as.factor(block)),
            position = position_dodge(width = 0.6)) +
  geom_point(aes(x = block, y = value, shape = name,
                 color = as.factor(block), group = as.factor(block)),
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
         shape = guide_legend(title.position="top", title.hjust = 0.5)) +
  coord_cartesian(ylim = c(0.5, 1))
  
ggsave(pl, file = "outputs/figures/FigS6-AUC.jpg", height = 5, width = 7)





# auc1 <- pivot_longer(auc, cols = !c("sp.code", "block", "source")) %>%
#   mutate(inout = case_when(name %in% c("AUCin", "AUCin.full", "in.cell", "in.full.cell", "in.full.n", "in.n") ~ "In sample",
#                            T ~ "Out of sample"),
#          type = case_when(name %in% c("AUCin", "AUCin.full", "AUCout", "AUCout.full") ~ "AUC",
#                           name %in% c("in.cell", "in.full.cell", "out.cell", "out.full.cell") ~ "cells",
#                           name %in% c("in.full.n", "in.n", "out.full.n", "out.n") ~ "samples"),
#          full = case_when(name %in% c("AUCin.full", "AUCout.full", "in.full.cell", "in.full.n", "out.full.cell", "out.full.n") ~ "full",
#                           T ~ "dataset")) %>%
#   select(!name) %>%
#   pivot_wider(names_from = c(type), values_from = value) %>%
#   mutate(source = case_when(full == "dataset" ~ source,
#                             T ~ "All"),
#          species = case_when(sp.code == "GPOR" ~ "Spring Salamander",
#                              sp.code == "RACA" ~ "Cascades Frog")) %>%
#   distinct()
# 
# pl <- ggplot(auc1) +
#   geom_line(aes(x = source, y = AUC,
#                 group = interaction(source, block), color = as.factor(block)),
#             position = position_dodge(width = 0.6)) +
#   geom_point(aes(x = source, y = AUC, shape = inout,
#                  color = as.factor(block), group = as.factor(block), size = cells),
#              position = position_dodge(width = 0.6)) +
#   theme_bw() +
#   scale_shape_manual(values = c(16, 1)) +
#   scale_color_manual(values = blockcols) +
#   labs(x = "Dataset", y = "AUC", color = "Excluded block",
#        shape = "Validation", size = "Number of cells \nwith validation data") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#         strip.background = element_blank(),
#         legend.position = "bottom") +
#   scale_x_discrete(labels = scales::label_wrap(15)) +
#   facet_wrap(~species, scales = "free_x") +
#   guides(color = guide_legend(title.position="top", title.hjust = 0.5),
#          size = guide_legend(title.position="top", title.hjust = 0.5),
#          shape = guide_legend(title.position="top", title.hjust = 0.5))
# 
# ggsave(pl, file = "outputs/figures/FigS6-AUC.jpg", height = 7, width = 11)





# Table S3: Data summary ----

dat <- c()

load("outputs/2_RACA_tau1/datafull-info.rdata")
dat1 <- species.data$locs$disc %>%
  mutate(species = "Cascades Frog",
         sp.code = "RACA")
dat <- bind_rows(dat, dat1)

load("outputs/8_GPOR_tau1/datafull-info.rdata")
dat1 <- species.data$locs$disc %>%
  mutate(species = "Spring Salamander",
         sp.code = "GPOR")
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
  arrange(species, data.type, desc(n.obs)) %>%
  distinct()
datsum$covariates[substr(datsum$covariates, 1, 1) == ","] <- gsub(", ", "", datsum$covariates[substr(datsum$covariates, 1, 1) == ","])
datsum$covariates[grep("^yday, $", datsum$covariates)] <- 'yday'

colnames(datsum) <- c("Species", "Data type", "Dataset", "Method", "Number of observations", 
                      "Number of sites", "Number of hexbins", "Extent", "Covariates")

write.csv(datsum, file = "outputs/tables/TabS3-datasummary.csv",
          row.names = F)




# Other: dataset intercepts ----


load("outputs/2_RACA_tau1/datafull.rdata")

dp <- out$alpha


a <- ggplot(dp) + 
  geom_pointrange(aes(x = reorder(name, mean), y = log(mean), 
                      ymax = log(hi), ymin = log(lo),
                      color = data.type)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1, 
                                   vjust = 0.5)) +
  labs(x = "Dataset", y = "Log(Dataset intercept)", color = "Data type") +
  coord_cartesian(ylim = c(-8, 3))




load("outputs/8_GPOR_tau1/datafull.rdata")

dp <- out$alpha


b <- ggplot(dp) + 
  geom_pointrange(aes(x = reorder(name, mean), y = log(mean), 
                      ymax = log(hi), ymin = log(lo),
                      color = data.type)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1, 
                                   vjust = 0.5),
        axis.text.y = element_blank()) +
  labs(x = "Dataset", y = "Log(Dataset intercept)", color = "Data type") +
  coord_cartesian(ylim = c(-8, 3))

pl <- a | b 
pl <- pl +
  plot_layout(guides = "collect",
              axes = "collect") +
  plot_annotation(tag_levels = "a") 
pl


# Appendix 3: tau ----

dirs <- list.dirs("outputs/", recursive = F)
dirs <- dirs[grep("RACA|GPOR", dirs)]

tau <- c()
auc <- c()
lambda <- c()
psi <- c()
for (d in 1:length(dirs)) {
  
  
  ind <- as.numeric(gsub("outputs/|_.*", "", dirs[d]))
  if (ind %in% c(1, 7)) mod <- 1
  if (ind %in% c(2, 8)) mod <- 2
  if (ind %in% c(3, 9)) mod <- 3
  if (ind %in% c(4, 10)) mod <- 4
  if (ind %in% c(5, 11)) mod <- 5
  if (ind %in% c(6, 12)) mod <- 6
  
  
  load(paste0(dirs[d], "/datafull-info.rdata"))
  
  auc1 <- all.auc %>%
    mutate(model = mod,
           block = "none")
  auc <- bind_rows(auc, auc1)
  
  for (i in 1:3) {
    load(paste0(dirs[d], "/data", i, "-info.rdata"))
    auc1 <- all.auc %>%
      mutate(model = mod,
             block = as.character(i))
    auc <- bind_rows(auc, auc1)
    
  }
  
  load(paste0(dirs[d], "/datafull.rdata"))
  
  tau1 <- out$tau %>%
    mutate(model = mod,
           sp.code = auc1$sp.code[1])
  
  tau <- bind_rows(tau, tau1)
  
  
  lambda1 <- out$lambda %>%
    mutate(model = mod,
           sp.code = auc1$sp.code[1])
  
  lambda <- bind_rows(lambda, lambda1)
  
  psi1 <- out$psi %>%
    mutate(model = mod,
           sp.code = auc1$sp.code[1])
  
  psi <- bind_rows(psi, psi1)
  
}


# Plot tau 
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
  facet_wrap(~sp.code + model, scales = "free_x") +
  coord_cartesian(ylim = c(-10, 10))

ggplot() + 
  geom_violin(data = priors, aes(x = name, y = prior)) +
  geom_pointrange(data = filter(tau, lo != hi), aes(x = name, y = mean, ymin = lo, ymax = hi)) +
  theme_bw() +
  labs(x = "Model", y = "Value") +
  facet_wrap(~sp.code + model, scales = "free_x") +
  coord_cartesian(ylim = c(-0.5, 0.5))



# Plot AUC
auc1 <- auc %>%
  select(sp.code, model, block, AUCin.full, AUCout.full) %>%
  distinct() %>%
  pivot_longer(cols = c(AUCin.full, AUCout.full))%>%
  mutate(dist = case_when(model == 4 ~ "Gamma(5, 5)",
                          model == 5 ~ "Normal(0, 3.16)",
                          model == 6 ~ "Gamma(0.01, 0.01)",
                          model == 1 ~ "0.1",
                          model == 2 ~ "1",
                          model == 3 ~ "10"),
         val = case_when(name == "AUCin.full" ~ "In-sample",
                         T ~ "Out-of-sample"))


blockcols <- c("none" = "black", "1" = "#e79f1e", "2" = "#009e73", "3" = "#cb79a8")
ggplot(auc1) +
  geom_point(aes(x = dist, y = value, shape = val, group = block, color = block),
             position = position_dodge(width = 0.6)) +
  geom_line(aes(x = dist, y = value, group = interaction(block, dist), color = block),
            position = position_dodge(width = 0.6)) +
  theme_bw() +
  scale_shape_manual(values = c(16, 1)) +
  scale_color_manual(values = blockcols) +
  labs(x = "Tau", y = "AUC", color = "Fold", shape = "Validation") +
  facet_wrap(~ sp.code)


# Maps
load("outputs/1_RACA_tau0.1/region.rdata")
lam <- filter(lambda, sp.code == "RACA") %>%
  full_join(region$sp.grid, ., by = "conus.grid.id")


q99 <- stats::quantile(lam$mean, 0.99, na.rm = T)
lam$mean[which(lam$mean > q99)] <- q99

ggplot(lam) +
  geom_sf(aes(fill = mean), color = NA) +
  facet_wrap(~ model) +
  viridis::scale_fill_viridis(guide = ggplot2::guide_colorbar(),
                              option = "magma",
                              na.value = NA,
                              direction = -1)




ggplot(lam) +
  geom_sf(aes(fill = log(mean)), color = NA) +
  facet_wrap(~ model) +
  viridis::scale_fill_viridis(guide = ggplot2::guide_colorbar(),
                              option = "magma",
                              na.value = NA,
                              direction = -1)


q99 <- stats::quantile(lam$mean, 0.95, na.rm = T)
lam$mean[which(lam$mean > q99)] <- q99




ggplot(lam) +
  geom_sf(aes(fill = mean), color = NA) +
  facet_wrap(~ model) +
  viridis::scale_fill_viridis(guide = ggplot2::guide_colorbar(),
                              option = "magma",
                              na.value = NA,
                              direction = -1)


tmp <- lam %>%
  st_drop_geometry() %>%
  select(conus.grid.id, model, mean) %>%
  pivot_wider(values_from = mean, names_from = model)
tmp1 <- as.data.frame(cor(tmp[,2:7])) %>%
  dplyr::mutate(model1 = row.names(.)) %>%
  tidyr::pivot_longer(!model1) %>%
  dplyr::filter(value != 1) 

ggplot(tmp1) + geom_tile(aes(x = model1, y = name, fill = value))

ps <- filter(psi, sp.code == "RACA") %>%
  mutate(pres = case_when(mean > 0.5 ~ "present",
                          T ~ "absent")) %>%
  full_join(region$sp.grid, ., by = "conus.grid.id")



ggplot(ps) +
  geom_sf(aes(fill = pres), color = NA) +
  facet_wrap(~ model) +
  ggplot2::scale_fill_manual(values = c("white", "darkblue"), na.value = NA) +
  ggplot2::scale_color_manual(values = c("white", "darkblue"), na.value = NA)



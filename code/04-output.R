

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



# Figure 1: Model framework ----



## Figure 2: RACA- range (a), intensity (b), suitable habitat (c), uncertainty (d) ----
# load(paste0(out.dir, "datafull.rdata"))
out.dir <- "outputs/1_RACA_iSDM/"
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

racamap <- maprange | mapocc | mapint | mapunc
racamap <- racamap + 
  plot_annotation(tag_levels = c("a"))
ggsave(racamap, file = "outputs/figures/Fig2-racamap.jpg",
       height = 7, width = 12)





## Figure 3: RACA- parameters ----
load("outputs/1_RACA_iSDM/datafull-info.rdata")
load("outputs/1_RACA_iSDM/datafull.rdata")
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
    
    # now get unscaled values
    mn <- mean(covar_unscaled[,cov])
    sd <- sd(covar_unscaled[,cov])
    
    use$x_unscaled <- (use$x * sd) + mn
    
    
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
    
    # now get unscaled values
    mn <- mean(covar_unscaled[,cov])
    sd <- sd(covar_unscaled[,cov])
    
    use$x_unscaled <- (use$x * sd) + mn
    
    
    all <- bind_rows(all, use)
  }
  

}


all1 <- inner_join(all, cov.labs, by = c("cov" = "covariate")) %>%
  mutate(cov = Label)


effects <- ggplot(filter(all1)) +
  geom_hline(yintercept = 0) +
  geom_ribbon(aes(x = x_unscaled, ymax = hi, ymin = lo), alpha = 0.5) +
  geom_line(aes(x = x_unscaled, y = mean)) +
  facet_wrap(~ cov, scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank())+
  labs(x = "Covariate value", y = "Exp(Estimate)")
effects


covs <- pars / effects
covs <- covs + 
  plot_annotation(tag_levels = "a") + 
  plot_layout(widths = c(1, 1.2))

ggsave(covs, file = "outputs/figures/Fig3-racaparameters.jpg",
       height = 6, width = 8)












# Figure 4: GPOR- Process model ----

load("outputs/2_GPOR_iSDM/datafull-info.rdata")
load("outputs/2_GPOR_iSDM/datafull.rdata")
load("outputs/2_GPOR_iSDM/region.rdata")

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
  geom_sf(data = region$range, fill = NA, color = "gray35", linewidth = 0.4, aes(linetype = range.name)) +
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
                               position = "bottom"),
         linetype = "none")


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


## Figure 5: GPOR- PO data and effort ----

load("outputs/2_GPOR_iSDM/datafull-info.rdata")
load("outputs/2_GPOR_iSDM/datafull.rdata")
load("outputs/2_GPOR_iSDM/region.rdata")




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

allfiles <- read.csv("data/00-data-summary-flexiSDM.csv") %>%
  filter(Species == "GPOR") %>%
  rename(file.name = Data.Swamp.file.name,
         file.label = Name,
         covar.mean = Covar.mean,
         covar.sum = Covar.sum,
         data.type = Type.true) %>%
  select(file.name, file.label, covar.mean, covar.sum, data.type, PO.extent)


species.data <- load_species_data("GPOR",
                                  sp.code.all,
                                  file.info = allfiles,
                                  file.path = "data/data-ready/",
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
  filter(source %in% c("VT PO", "NY Atlas", "WV PO", "MA PO", "MD PO", "ME PO")) %>%
  mutate(map = "state")
use1$source <- factor(use1$source, levels = c("VT PO", "NY Atlas", "WV PO", "MA PO", "MD PO", "ME PO", "iNaturalist"))

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
                                "MD PO" = "#E63833FF",
                                "ME PO" = "#881C00FF",
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
use2$source <- factor(use2$source, levels = c("VT PO", "NY Atlas", "WV PO", "MA PO", "MD PO", "ME PO", "iNaturalist"))

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
                                "MD PO" = "#E63833FF",
                                "ME PO" = "#881C00FF",
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
  filter(PO.dataset.name %in% c("MA PO, MD PO, ME PO, NY Atlas, VT PO, WV PO")) %>%
  full_join(region$sp.grid, ., by = "conus.grid.id") %>%
  filter(mean != 0)

pl3 <- base +
  geom_sf(data = tmp, aes(fill = log(mean), color = log(mean))) +
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


tmp1 <- out$effort %>%
  filter(PO.dataset.name %in% c("iNaturalist")) %>%
  full_join(region$sp.grid, ., by = "conus.grid.id") 
tmp <- tmp1 %>%
  filter(mean != 0)

pl4 <- base +
  geom_sf(data = tmp1, fill = "white", color = 'white') +
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
  filter(name %in% c("MA PO, MD PO, ME PO, NY Atlas, VT PO, WV PO")) %>%
  mutate(covariate = case_when(covariate == "traveltime" ~ "Travel time",
                               T ~ covariate))
pars1$covariate <- factor(pars1$covariate, levels = c("MD", "ME", "NY", "VT", "WV", "Travel time"))

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
  mutate(covariate = case_when(covariate == "ORM" ~ "ORM",
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

ggsave(pl, file = "outputs/figures/Fig5-gporPO.jpg",
       height = 12, width = 12)








# Figure S2: CV maps ----

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
load("outputs/1_RACA_iSDM/datafull-info.rdata")
load("outputs/1_RACA_iSDM/region.rdata")

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
load("outputs/2_GPOR_iSDM/datafull-info.rdata")
load("outputs/2_GPOR_iSDM/region.rdata")

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
load("outputs/2_GPOR_iSDM/region.rdata")
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


out.dir <- "outputs/1_RACA_iSDM/"
# Load data from all model runs
for (i in 1:3) {
  load(paste0(out.dir, "data", i, "-info.rdata"))
  all.auc$block <- as.character(all.auc$block)
  auc <- bind_rows(auc, all.auc)
}
# now add full model
load(paste0(out.dir, "datafull-info.rdata"))
auc <- bind_rows(auc, all.auc)

out.dir <- "outputs/2_GPOR_iSDM/"
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
  mutate(species = case_when(sp.code == "GPOR" ~ "Spring salamander",
                             sp.code == "RACA" ~ "Cascades frog"))

mean(auc1$AUCout.full, na.rm = T)
mean(auc1$AUCin.full, na.rm = T)

auc2 <- pivot_longer(auc1, cols = !c("species", "sp.code", "block")) %>%
  mutate(inout = case_when(name == "AUCin.full" ~ "In-sample",
                           name == "AUCout.full" ~ "Out-of-sample"))

pla <- ggplot(filter(auc2, species == "Cascades frog")) +
  geom_hline(yintercept = 0.5, color = "darkgray") +
  geom_line(aes(x = block, y = value,
                group = interaction(block), color = as.factor(block)),
            position = position_dodge(width = 0.6)) +
  geom_point(aes(x = block, y = value, shape = inout,
                 color = as.factor(block), group = as.factor(block)),
             position = position_dodge(width = 0.6)) +
  theme_bw() +
  scale_shape_manual(values = c(16, 1)) +
  scale_color_manual(values = blockcols) +
  labs(x = "Dataset", y = "AUC", color = "Excluded fold",
       shape = "Validation") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        strip.background = element_blank(),
        legend.position = "bottom") +
  scale_x_discrete(labels = scales::label_wrap(15)) +
  #facet_wrap(~species, scales = "free_x") +
  guides(color = guide_legend(title.position="top", title.hjust = 0.5),
         shape = guide_legend(title.position="top", title.hjust = 0.5)) +
  coord_cartesian(ylim = c(0.4, 1))

plb <- ggplot(filter(auc2, species == "Spring salamander")) +
  geom_hline(yintercept = 0.5, color = "darkgray") +
  geom_line(aes(x = block, y = value,
                group = interaction(block), color = as.factor(block)),
            position = position_dodge(width = 0.6)) +
  geom_point(aes(x = block, y = value, shape = inout,
                 color = as.factor(block), group = as.factor(block)),
             position = position_dodge(width = 0.6)) +
  theme_bw() +
  scale_shape_manual(values = c(16, 1)) +
  scale_color_manual(values = blockcols) +
  labs(x = "Dataset", y = "AUC", color = "Excluded fold",
       shape = "Validation") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        strip.background = element_blank(),
        legend.position = "bottom") +
  scale_x_discrete(labels = scales::label_wrap(15)) +
  #facet_wrap(~species, scales = "free_x") +
  guides(color = guide_legend(title.position="top", title.hjust = 0.5),
         shape = guide_legend(title.position="top", title.hjust = 0.5)) +
  coord_cartesian(ylim = c(0.4, 1))
  


pl <- pla | plb 
pl <- pl +
  plot_layout(guides = "collect",
              axes = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom")

pl
ggsave(pl, file = "outputs/figures/FigS6-AUC.jpg", height = 4, width = 7)








# Figure S7: Dataset intercepts ----


load("outputs/1_RACA_iSDM/datafull.rdata")

dp <- out$alpha %>%
  arrange(mean)
dp$name <- factor(dp$name, levels = dp$name)

a <- ggplot(dp) + 
  geom_pointrange(aes(x = name, y = log(mean), 
                      ymax = log(hi), ymin = log(lo),
                      color = data.type)) +
  scale_color_manual(values = c("#000000", "#56B4E9", "#CC79A7")) +
  scale_x_discrete(labels = stringr::str_wrap(dp$name, 20)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1, 
                                   vjust = 0.5),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "Dataset", y = "Log(Dataset intercept)", color = "Data type") +
  coord_cartesian(ylim = c(-10, 3))




load("outputs/2_GPOR_iSDM/datafull.rdata")

dp <- out$alpha %>%
  arrange(mean)
dp$name <- factor(dp$name, levels = dp$name)


b <- ggplot(dp) + 
  geom_pointrange(aes(x = name, y = log(mean), 
                      ymax = log(hi), ymin = log(lo),
                      color = data.type)) +
  scale_color_manual(values = c("#000000", "#56B4E9", "#CC79A7")) +
  scale_x_discrete(labels = stringr::str_wrap(dp$name, 20)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1, 
                                   vjust = 0.5),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "Dataset", y = "Log(Dataset intercept)", color = "Data type") +
  coord_cartesian(ylim = c(-10, 3))

b
pl <- a | b 
pl <- pl +
  plot_layout(guides = "collect",
              axes = "collect") +
  plot_annotation(tag_levels = "a") 

ggsave(pl, file = "outputs/figures/FigS7-datasetintercepts.jpg", height = 5, width = 11)




# Figure S8: GPOR marginal effects ----

load("outputs/2_GPOR_iSDM/datafull-info.rdata")
load("outputs/2_GPOR_iSDM/datafull.rdata")
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
    
    # now get unscaled values
    mn <- mean(covar_unscaled[,cov])
    sd <- sd(covar_unscaled[,cov])
    
    use$x_unscaled <- (use$x * sd) + mn
    
    
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
    
    # now get unscaled values
    mn <- mean(covar_unscaled[,cov])
    sd <- sd(covar_unscaled[,cov])
    
    use$x_unscaled <- (use$x * sd) + mn
    
    
    all <- bind_rows(all, use)
  }
  
  
}


all1 <- inner_join(all, cov.labs, by = c("cov" = "covariate")) %>%
  mutate(cov = Label,
         x_unscaled = case_when(cov == "Forest (% cover)" ~ x_unscaled * 100,
                          T ~ x_unscaled))


effects <- ggplot(filter(all1)) +
  geom_hline(yintercept = 0) +
  geom_ribbon(aes(x = x_unscaled, ymax = hi, ymin = lo), alpha = 0.5) +
  geom_line(aes(x = x_unscaled, y = mean)) +
  facet_wrap(~ cov, scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank())+
  labs(x = "Covariate value", y = "Exp(Estimate)")


effects
ggsave(effects, file = "outputs/figures/FigS8-GPORmarginaleffects.jpg",
       height = 7, width = 7)






# Table S3: Data summary ----

dat <- c()

load("outputs/1_RACA_iSDM/datafull-info.rdata")
dat1 <- species.data$locs$disc %>%
  mutate(species = "Cascades Frog",
         sp.code = "RACA")
dat <- bind_rows(dat, dat1)

load("outputs/2_GPOR_iSDM/datafull-info.rdata")
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



